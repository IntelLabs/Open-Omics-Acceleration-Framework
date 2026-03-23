# Multiprocess Task Runner

The Multiprocess Task Runner is a generic utility designed to improve CPU utilization by launching multiple parallel processes for workloads that do not scale efficiently on a single process. It reads a JSON configuration file describing your tasks and executes them in parallel with automatically managed multiprocessing.

The tool supports multiple common execution patterns and manages parallelism automatically.

---

## General Usage

Run the script using the following command:

```bash
python multiprocess.py --json_file=<configuration.json> --case=<case_number>
```

### Optional Parameters

- **`--single_socket`**: Restricts execution to a single socket.
- **`--include_hyperthreads`**: Enables hyperthreads during execution.

### Viewing Help

For more details about the available options, run the following command:

```bash
python multiprocess.py --help
```

## Case 1: Tasks Defined by Input Files

For this case, the tasks are determined by input files. Your configuration JSON file should include:

```json
{
    "script_command": "python myTool.py",
    "input_dir": "my_input/",
    "output_log_dir": "output_log/",
    "shared_args": "disable_cuda=True noipex=True bf16=False timing=True",
    "input_file_arg": "input_file="
}
```

To execute, run:

```bash
python multiprocess.py --json_file=config.json --case=1
```

## Case 2: Tasks Defined by Unique Arguments

Here, the tasks are determined by a set of unique arguments. Your configuration JSON file should include:

```json
{
    "script_command": "python myTool.py",
    "output_log_dir": "output_log/",
    "shared_args": "disable_cuda=True timing=True",
    "unique_args": [
        "num_iter=100 noipex=True bf16=False",
        "num_iter=200 noipex=False bf16=True",
        "num_iter=300 noipex=False bf16=True"
    ]
}
```

To execute, run:

```bash
python multiprocess.py --json_file=config.json --case=2
```

## Case 3: Tasks Defined by Input Files and Unique Arguments

For this case, each task is defined by a combination of input files and unique arguments. Your configuration JSON file should include:

```json
{
    "script_command": "python myTool.py",
    "input_dir": "my_input/",
    "output_log_dir": "output_log/",
    "shared_args": "disable_cuda=True timing=True",
    "input_file_arg": "input_file=",
    "tasks": [
        {
            "input_file": "file1.txt",
            "unique_args": "num_iter=100 noipex=True bf16=False"
        },
        {
            "input_file": "file2.txt",
            "unique_args": "num_iter=200 noipex=False bf16=True"
        }
    ]
}
```

To execute, run:

```bash
python multiprocess.py --json_file=config.json --case=3
```

## Additional Settings for Case 1 and Case 3

For Case 1 and Case 3, you can optionally include the following settings in your configuration JSON:

- sort: A boolean value (true/false) to specify if tasks should be processed in sorted order of input file size.
- reverse: A boolean value (true/false) to specify the sorting order:
  - true: Descending order
  - false: Ascending order

By default, tasks are processed in descending order of file size.
