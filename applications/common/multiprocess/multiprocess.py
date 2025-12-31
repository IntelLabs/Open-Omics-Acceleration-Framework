import re
import subprocess
import os
import psutil
import time
import multiprocessing as mp
from datetime import datetime
from pathlib import Path
import json
from absl import app, flags

# Define flags
FLAGS = flags.FLAGS

# Positional arguments
flags.DEFINE_string(
    "json_file", None, "Path to the JSON file containing task configurations."
)
flags.DEFINE_integer(
    "case", None, "Case type: 1, 2, or 3.", lower_bound=1, upper_bound=3
)

# Optional arguments
flags.DEFINE_bool("single_socket", False, "Use only a single socket. Default is false.")
flags.DEFINE_bool("include_hyperthreads", False, "Use hyperthreads. Default is false.")

# Mark positional arguments as required
flags.mark_flag_as_required("json_file")
flags.mark_flag_as_required("case")


class Task:

    used_identifiers = set()

    def __init__(self, command, unique_id, output_dir):
        if unique_id in Task.used_identifiers:
            raise ValueError(f"Task identifier '{unique_id}' is not unique.")

        self.command = command
        self.unique_id = unique_id
        self.output_log_dir = output_dir

        # Add the task identifier to the used_identifiers set
        Task.used_identifiers.add(unique_id)

    def __repr__(self):
        return f"Task(command={self.command}, identifier={self.unique_id}, output_dir={self.output_log_dir})"


class SystemInfo:
    def __init__(self, single_socket=False, include_hyperthreads=False):
        """
        Initialize the system information.

        Args:
            single_socket (bool): If True, gather information for a single socket only.
            include_hyperthreads (bool): If True, include hyperthreads in calculations.
        """
        self.single_socket = single_socket
        self.include_hyperthreads = include_hyperthreads
        self.num_sockets = self._get_num_sockets()
        self.total_numa_nodes, self.numa_nodes_list = self._get_numa_nodes()
        self.total_cores, self.core_list = self._get_core_list()

    def _get_num_sockets(self):
        """
        Determine the number of sockets on the current machine.

        Returns:
            int: Number of sockets.
        """
        try:
            output = subprocess.check_output("lscpu", shell=True, text=True)
            for line in output.splitlines():
                if "Socket(s):" in line:
                    return int(line.split(":")[1].strip())
        except Exception as e:
            raise RuntimeError("Unable to determine the number of sockets.") from e

    def _get_numa_nodes(self):
        """
        Determine the number of NUMA nodes, adjusted for single-socket configuration.

        Returns:
            list: A list of NUMA nodes to consider.
        """
        try:
            output = subprocess.check_output("lscpu", shell=True, text=True)
            total_numa_nodes = 0
            for line in output.splitlines():
                if "NUMA node(s):" in line:
                    total_numa_nodes = int(line.split(":")[1].strip())

            if self.single_socket:
                # If we are considering single socket, divide NUMA nodes equally across all sockets
                nodes_per_socket = total_numa_nodes // self.num_sockets
                numa_nodes = list(range(nodes_per_socket))
            else:
                numa_nodes = list(range(total_numa_nodes))

            return len(numa_nodes), numa_nodes
        except Exception as e:
            raise RuntimeError("Unable to determine the number of NUMA nodes.") from e

    def _get_core_list(self):
        """
        Create a core list and calculate the total cores, adjusted for hyperthreading and NUMA distribution.

        Returns:
            tuple: Total cores and core list.
        """
        try:
            # Get the lscpu output with NUMA node CPU info
            output = subprocess.check_output("lscpu", shell=True, text=True)

            numa_cores = {}
            # Regular expression to find NUMA nodes and their corresponding CPUs
            numa_pattern = re.compile(r"NUMA node(\d+)\s+CPU\(s\):\s+([\d,-]+)")

            for line in output.splitlines():
                match = numa_pattern.search(line)
                if match:
                    node_id = int(match.group(1))  # NUMA node ID
                    cpu_range = match.group(2)  # CPU range for this NUMA node

                    # Process the CPU range for each node
                    cores = []
                    cpu_range = cpu_range.split(",")
                    if not self.include_hyperthreads:
                        cpu_range = [cpu_range[0]]

                    for cpu_range_part in cpu_range:
                        start, end = map(int, cpu_range_part.split("-"))
                        cores.extend(range(start, end + 1))
                    numa_cores[node_id] = cores

            # Filter NUMA nodes based on the 'single_socket' flag
            if self.single_socket:
                # Get the total number of CPUs available for all NUMA nodes within a socket
                # Assume that each socket has equal NUMA nodes
                socket_size = len(numa_cores) // self.num_sockets
                nodes_to_consider = list(numa_cores.keys())[:socket_size]
            else:
                nodes_to_consider = list(numa_cores.keys())

            # Now, ensure that each NUMA node contributes an equal number of cores
            num_cores_per_node = min(len(cores) for cores in numa_cores.values())
            core_list = []
            for node_id in nodes_to_consider:
                # Take the first `num_cores_per_node` cores from each NUMA node
                core_list.extend(numa_cores[node_id][:num_cores_per_node])

            total_cores = len(core_list)
            return total_cores, core_list
        except Exception as e:
            raise RuntimeError("Unable to create the core list.") from e

    def __str__(self):
        """
        String representation of the system information.

        Returns:
            str: A formatted string with system information.
        """
        return (
            f"System Information:\n"
            f"Single Socket: {self.single_socket}\n"
            f"Include Hyperthreads: {self.include_hyperthreads}\n"
            f"Total NUMA Nodes: {self.total_numa_nodes}\n"
            f"NUMA Nodes: {self.numa_nodes_list}\n"
            f"Total Cores: {self.total_cores}\n"
            f"Core List: {self.core_list}"
        )


def output_dir_create(output_dir):
    # Generate a unique directory name with a timestamp
    unique_dir_name = f"{datetime.now().strftime('%Y_%m_%d_%H%M%S')}_"

    # Get the current working directory and construct the full path
    curr_output_dir  = os.path.join(os.getcwd(), output_dir, unique_dir_name)

    # Create the directory (and intermediate directories if they don't exist)
    os.makedirs(curr_output_dir , exist_ok=True)

    print(f"Output log directory created: {curr_output_dir}")

    return curr_output_dir

def start_bash_subprocess(task: Task, mem, core_list):
    """Starts a new bash subprocess and puts it on the specified cores."""

    os.makedirs(task.output_log_dir, exist_ok=True)
    log_file_name = os.path.join(task.output_log_dir, str(task.unique_id) + ".txt")

    base_fold_cmd = "/usr/bin/time -v {} "
    command = base_fold_cmd.format(task.command)
    numactl_args = [
        "numactl",
        "-m",
        mem,
        "-C",
        "-".join([str(core_list[0]), str(core_list[-1])]),
        command,
    ]
    command = " ".join(numactl_args)
    print(command)

    try:
        with open(log_file_name, "a") as f:
            process = subprocess.call(
                command, shell=True, universal_newlines=True, stdout=f, stderr=f
            )
    except Exception as e:
        print("exception for unique task=", str(task.unique_id), e)
        process = None
    return (process, mem, task, core_list)


def check_available_memory():
    """Checks for available memory using psutil."""
    mem = psutil.virtual_memory()
    available_memory = mem.available
    return available_memory / 1024**2


def multiprocessing_run(tasks, max_processes, system_info: SystemInfo):
    core_list = system_info.core_list
    cores_per_process = system_info.total_cores // max_processes
    pool = mp.Pool(processes=max_processes)

    queue = [i for i in range(max_processes)]
    error_tasks = []

    def update_queue(result):
        queue.append(result[3][0] // cores_per_process)
        if result[0] != 0:
            error_tasks.append(result[2])

    # Iterate over the files and start a new subprocess for each file.
    print("Total Tasks = ", len(tasks))
    results = [None] * len(tasks)

    # numa_nodes
    numa_nodes = system_info.total_numa_nodes

    i = 0
    for task in tasks:

        process_num = queue.pop(0)

        if max_processes == 1:
            if numa_nodes > 1:
                mem = "0-{}".format(numa_nodes - 1)
            else:
                mem = "0"
        else:
            mem = str(process_num // (max_processes // numa_nodes))

        results[i] = pool.apply_async(
            start_bash_subprocess,
            args=(
                task,
                mem,
                core_list[
                    process_num
                    * cores_per_process : (process_num + 1)
                    * cores_per_process
                ],
            ),
            callback=update_queue,
        )
        i += 1
        while len(queue) == 0 and i < len(tasks):
            time.sleep(0.05)

    pool.close()
    pool.join()

    return error_tasks


def combine_command(
    script_command: str,
    shared_args: str,
    unique_args: str = "",
    input_file_arg: str = "",
    input_file: str = "",
) -> str:
    """
    Combines script command, shared arguments, unique arguments, and input file details into a complete command string.
    Supports various combinations of arguments based on different cases.

    Args:
        script_command (str): The base script command.
        shared_args (str): Shared arguments for the script.
        unique_args (str, optional): Unique arguments for the script. Defaults to "".
        input_file_arg (str, optional): The input file argument prefix. Defaults to "".
        input_file (str, optional): The input file. Defaults to "".

    Returns:
        str: The complete command string.
    """
    import shlex

    # Split components using shlex to handle any spaces or special characters in the arguments
    script_command_list = shlex.split(script_command)
    shared_args_list = shlex.split(shared_args)
    unique_args_list = shlex.split(unique_args)

    # Start building the command with the script command and shared arguments
    command_parts = script_command_list + shared_args_list

    # Add unique arguments if provided
    if unique_args:
        command_parts += unique_args_list

    # Add input file argument and file if provided
    if input_file_arg and input_file:
        command_parts.append(f"{input_file_arg}{input_file}")

    # Join all parts into a final command string
    final_command = shlex.join(command_parts)

    return final_command


def validate_json(config, case):
    """Validate the JSON configuration based on the case."""
    required_keys = []

    if case == 1:
        required_keys = [
            "script_command",
            "input_dir",
            "output_log_dir",
            "shared_args",
            "input_file_arg",
        ]
    elif case == 2:
        required_keys = [
            "script_command",
            "output_log_dir",
            "shared_args",
            "unique_args",
        ]
    elif case == 3:
        required_keys = [
            "script_command",
            "input_dir",
            "output_log_dir",
            "shared_args",
            "input_file_arg",
            "tasks",
        ]

    missing_keys = [key for key in required_keys if key not in config]
    if missing_keys:
        raise ValueError(
            f"Missing required keys for case {case}: {', '.join(missing_keys)}"
        )

    # Optional keys with default values
    if case in [1, 3]:
        config.setdefault("sort", True)
        config.setdefault("reverse", True)


def common_preparation(config):
    output_dir = config["output_log_dir"]
    output_log_dir = output_dir_create(output_dir)
    config["output_log_dir"] = output_log_dir


def get_file_size(file_path):
    """Gets the size of the file in bytes."""
    return os.path.getsize(file_path)


def sort_files(input_files, reverse=True):
    """Sorts files based on their sizes."""
    size_dict = {file: get_file_size(file) for file in input_files}
    sorted_files = sorted(size_dict.keys(), key=lambda x: size_dict[x], reverse=reverse)
    return sorted_files


def get_input_files(input_dir, sort=True, reverse=True):

    # Validate input directory
    if not os.path.isdir(input_dir):
        raise ValueError(f"Input directory does not exist: {input_dir}")

    """Get sorted files from a directory."""
    input_files = [
        os.path.abspath(os.path.join(input_dir, f))
        for f in os.listdir(input_dir)
        if os.path.isfile(os.path.join(input_dir, f))
    ]

    if sort:
        input_files = sort_files(input_files, reverse=reverse)

    return input_files


def prepare_tasks_case1(config):
    """Prepare tasks for case 1: Tasks defined by input files."""
    tasks = []
    input_dir = config["input_dir"]
    output_dir = config["output_log_dir"]
    script_command = config["script_command"]
    shared_args = config["shared_args"]
    input_file_arg = config["input_file_arg"]

    sort = config.get("sort", True)
    reverse = config.get("reverse", True)

    input_files = get_input_files(input_dir, sort=sort, reverse=reverse)

    for input_file in input_files:
        curr_command = combine_command(
            script_command=script_command,
            shared_args=shared_args,
            input_file_arg=input_file_arg,
            input_file=input_file,
        )
        t = Task(
            command=curr_command, unique_id=Path(input_file).stem, output_dir=output_dir
        )
        tasks.append(t)

    return tasks


def prepare_tasks_case2(config):
    """Prepare tasks for case 2: Tasks defined by argument sets."""
    tasks = []
    output_dir = config["output_log_dir"]
    script_command = config["script_command"]
    shared_args = config["shared_args"]

    task_id = 0
    for unique_args in config["unique_args"]:
        curr_command = combine_command(
            script_command=script_command,
            shared_args=shared_args,
            unique_args=unique_args,
        )
        t = Task(command=curr_command, unique_id=str(task_id), output_dir=output_dir)

        print("Adding task =", t)
        tasks.append(t)
        task_id += 1
    return tasks


def prepare_tasks_case3(config):
    """Prepare tasks for case 3: Tasks defined by input files and argument sets."""
    tasks = []
    input_dir = config["input_dir"]
    output_dir = config["output_log_dir"]
    script_command = config["script_command"]
    shared_args = config["shared_args"]
    input_file_arg = config["input_file_arg"]

    sort = config.get("sort", True)
    reverse = config.get("reverse", True)

    # Validate input directory
    if not os.path.isdir(input_dir):
        raise ValueError(f"Input directory does not exist: {input_dir}")

    # todo: sorting not implemented

    task_id = 0
    for task in config["tasks"]:
        input_file = os.path.join(input_dir, task["input_file"])

        # Validate input file
        if not os.path.isfile(input_file):
            raise ValueError(f"Input file does not exist: {input_file}")

        unique_args = task.get("unique_args", "")

        curr_command = combine_command(
            script_command=script_command,
            shared_args=shared_args,
            input_file_arg=input_file_arg,
            input_file=input_file,
            unique_args=unique_args,
        )
        t = Task(command=curr_command, unique_id=str(task_id), output_dir=output_dir)

        print("Adding task =", t)
        tasks.append(t)
        task_id += 1

    return tasks


def main(argv):

    t1 = time.time()

    # Initialize System Info
    system_info = SystemInfo(
        single_socket=FLAGS.single_socket,
        include_hyperthreads=FLAGS.include_hyperthreads,
    )
    print(system_info)

    # Check if the JSON file exists
    if not os.path.isfile(FLAGS.json_file):
        raise FileNotFoundError(f"JSON file does not exist: {FLAGS.json_file}")

    # Load configurations from the JSON file
    with open(FLAGS.json_file, "r") as f:
        config = json.load(f)

    # Validate JSON structure
    validate_json(config, FLAGS.case)

    common_preparation(config)

    # Prepare tasks based on the case
    if FLAGS.case == 1:
        tasks = prepare_tasks_case1(config)
    elif FLAGS.case == 2:
        tasks = prepare_tasks_case2(config)
    elif FLAGS.case == 3:
        tasks = prepare_tasks_case3(config)

    # Print JSON file content in organized manner
    print("Loaded JSON Configuration:")
    print(json.dumps(config, indent=4))

    """The main function."""

    # numa_nodes
    cores_per_numa = system_info.total_cores // system_info.total_numa_nodes

    if cores_per_numa % 48 == 0:
        max_processes_list = [
            48 * system_info.total_numa_nodes,
            24 * system_info.total_numa_nodes,
            12 * system_info.total_numa_nodes,
            6 * system_info.total_numa_nodes,
            3 * system_info.total_numa_nodes,
            system_info.total_numa_nodes,
            1,
        ]
    elif cores_per_numa % 8 == 0:
        max_processes_list = [
            8 * system_info.total_numa_nodes,
            4 * system_info.total_numa_nodes,
            2 * system_info.total_numa_nodes,
            system_info.total_numa_nodes,
            1,
        ]
    elif cores_per_numa % 4 == 0:
        max_processes_list = [
            4 * system_info.total_numa_nodes,
            2 * system_info.total_numa_nodes,
            system_info.total_numa_nodes,
            1,
        ]
    elif cores_per_numa % 3 == 0:
        max_processes_list = [
            3 * system_info.total_numa_nodes,
            system_info.total_numa_nodes,
            1,
        ]
    elif cores_per_numa % 2 == 0:
        max_processes_list = [
            2 * system_info.total_numa_nodes,
            system_info.total_numa_nodes,
            1,
        ]
    else:
        max_processes_list = [system_info.total_numa_nodes, 1]

    print("Total cores: ", system_info.total_cores)
    print("Total memory: {} MB ".format(check_available_memory()))

    for max_processes in max_processes_list:
        os.environ["OMP_NUM_THREADS"] = str(system_info.total_cores // max_processes)
        print(
            "Number of OMP Threads = {}, for {} instances".format(
                os.environ.get("OMP_NUM_THREADS"), max_processes
            )
        )
        if len(tasks) >= max_processes:
            returned_tasks = multiprocessing_run(tasks, max_processes, system_info)
            if len(returned_tasks) > 0:
                print(
                    "Following task couldn't be processed with {} instances".format(
                        max_processes
                    )
                )
                for failed_task in returned_tasks:
                    print(failed_task)
            tasks = returned_tasks
        else:
            continue

    if len(tasks) > 0:
        print("Following tasks couldn't be processed")
        for task in tasks:
            print(task)

    t2 = time.time()
    print("Total Time in seconds = ", (t2 - t1))


if __name__ == "__main__":
    app.run(main)
