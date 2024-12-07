import json
import subprocess
from subprocess import Popen, PIPE, run
import os, sys, time
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def HWConfigure(sso, num_nodes, th=20):
    run('lscpu > lscpu.txt', capture_output=True, shell=True)
    dt={}
    with open('lscpu.txt', 'r') as f:
        l = f.readline()
        while l:
            try:
                a,b = l.strip('\n').split(':')
                #aa, bb = a.split(' ')
                #print(a, b)
                dt[a] = b
            except:
                pass
            
            l = f.readline()

    ncpus = int(dt['CPU(s)'])
    nsocks = int(dt['Socket(s)'])
    nthreads = int(dt['Thread(s) per core'])
    ncores = int(dt['Core(s) per socket'])
    nnuma = int(dt['NUMA node(s)'])

    if sso:
        nsocks = 1

    #th = 16  ## max cores per rank
    #num_nodes = 1
        
    num_physical_cores_all_nodes = num_nodes * nsocks * ncores
    num_physical_cores_per_node = nsocks * ncores
    num_physical_cores_per_rank = nsocks * ncores
    
    while num_physical_cores_per_rank > th:
        num_physical_cores_per_rank /= 2

    num_physical_cores_per_rank = int(num_physical_cores_per_rank)
    assert num_physical_cores_per_rank > 8, 'cores per rank should be > 8'
    
    N = int(num_physical_cores_all_nodes / num_physical_cores_per_rank)
    PPN = int(num_physical_cores_per_node / num_physical_cores_per_rank)
    CPUS = int(ncores * nthreads * nsocks / PPN - 2*nthreads)
    THREADS = CPUS
    #print(f"N={int(N)}")
    #print(f"PPN={int(PPN)}")
    #print(f"CPUS={int(CPUS)}")
    #print(f"THREADS={int(THREADS)}")
    
    threads_per_rank = num_physical_cores_per_rank * nthreads
    bits = pow(2, num_physical_cores_per_rank) - 1
    allbits = 0
    mask="["
    for r in range(N):
        allbits = allbits | (bits << r*num_physical_cores_per_rank)
        allbits = allbits | (allbits << nsocks * ncores)
        #print("{:x}".format(allbits))
        if mask == "[":
            mask = mask + hex(allbits)
        else:
            mask = mask+","+ hex(allbits)
        allbits=0
        #print("{:x}".format(mask))
    mask=mask + "]"
    #print("I_MPI_PIN_DOMAIN={}".format(mask))

    return N, PPN, CPUS, THREADS, mask



if __name__ == '__main__':
    ## rgs parser
    parser=ArgumentParser()
    parser.add_argument('--ref', default="", help="Input data directory")
    parser.add_argument('--reads', nargs='+', help="reads, expects both the reads at the same location")
    #parser.add_argument('--input', default="/input", help="Input data directory")
    parser.add_argument('--tempdir',default="",help="Intermediate data directory")
    #parser.add_argument('--refdir',default="/refdir",help="Reference genome directory")
    parser.add_argument('--output',default="/output", help="Output data directory")
    #parser.add_argument("-i", "--refindex", default="None", help="name of refindex file")
    #parser.add_argument("-r1", "--read1", default="None",  help="name of read1")
    #parser.add_argument("-r2", "--read2", default="None",  help="name of read2")
    parser.add_argument('--read3',default="None",  \
                        help="name of r3 files (for fqprocess) seperated by spaces")
    parser.add_argument('--readi1',default="None", \
                        help="name of i1 files (for fqprocess) seperated by spaces")
    parser.add_argument('--prefix',default="None", \
                        help="prefix for processed R1 and R3 files for bwa-mem2")
    parser.add_argument('--suffix',default="None", \
                        help="suffix for processed R1 and R3 files for bwa-mem2")   
    parser.add_argument('--whitelist',default="whitelist.txt",help="10x whitelist file")
    parser.add_argument('--read_structure',default="16C",help="read structure")
    parser.add_argument('--barcode_orientation',default="FIRST_BP_RC",help="barcode orientation")
    parser.add_argument('--sample_id',default=-1,help="sample id")
    parser.add_argument('--output_format',default="None",help="output_format")
    parser.add_argument("-b", "--bam_size",default=9, help="bam file size in GB")
    parser.add_argument('--mode', default='sortedbam', \
                        help="flatmode/fqprocessonly/multifq/sortedbam. flatmode is just bwa w/o sort.")
    parser.add_argument('--read_type',default="short",
                        help="(short/long): bwa-mem2 with short reads, mm2-fast with long reads")
    parser.add_argument('-in', '--rindex',action='store_true',help="It will index reference genome for bwa-mem2. If it is already done offline then don't use this flag.")
    parser.add_argument('--dindex',action='store_true',help="It will create .fai index. If it is done offline then disable this.")
    #parser.add_argument('--container_tool',default="docker",help="Container tool used in pipeline : Docker/Podman")
    parser.add_argument('-pr', '--profile',action='store_true',help="Use profiling")
    parser.add_argument('--keep_unmapped',action='store_true',help="keep unmapped entries at the end of sam file.")
    parser.add_argument('--params', type=str, default='@RG\\tID:RG1\\tSM:RGSN1', help="parameter string to bwa-mem2 barring threads paramter")
    #parser.add_argument("-p", "--outfile", default="final", help="prefix for read files")
    parser.add_argument("--sso", action='store_true', help="prefix for read files")
    parser.add_argument("--th", default=20, help="#min cores per rank")
    parser.add_argument("-N", default=-1, help="#ranks")
    parser.add_argument("-PPN", default=-1, help="ppn")
    parser.add_argument("--cpus", default=-1, help="CPUS")
    #parser.add_argument("--threads", default=-1, help="THREADS")
    args = vars(parser.parse_args())

    assert len(args["reads"]) >= 1
    
    args["input"] = os.path.dirname(args["reads"][0])

    args["read1"] = os.path.basename(args["reads"][0])
    if len(args["reads"]) == 2: args["read2"] = os.path.basename(args["reads"][1])
    else: args["read2"] = ""

    args["refdir"] = os.path.dirname(args["ref"])
    args["refindex"] = os.path.basename(args["ref"])
    args["output"] = os.path.dirname(args["output"])
    args["outfile"] = os.path.basename(args["output"])
    
    num_nodes=1
    N, PPN, CPUS, THREADS, mask = HWConfigure(args["sso"], num_nodes, args['th'])
    if args["N"] != -1: N = args["N"]
    if args["PPN"] != -1: PPN = args["PPN"]
    if args["cpus"] != -1: CPUS = args["cpus"]
    #if args["threads"] != -1: THREADS = args["cpus"]
    
    print("[Info] Running {} processes per compute node, each with {} threads".format(N, THREADS))
    args['cpus'], args['threads'] = str(CPUS), str(THREADS)

    cmd="hostname > hostfile"
    a = run(cmd, capture_output=True, shell=True)
    
    BINDING="socket"
    cmd="mkdir -p logs"
    a = run(cmd, capture_output=True, shell=True)

    #cmd = "export I_MPI_PIN_DOMAIN==mask" + "; mpiexec -bootstrap ssh -n " + N + "-ppn " + PPN + " -bind-to " + BINDING + "-map-by " + BINDING + " --hostfile hostfile  python -u fq2bams.py --cpus" + CPUS + " --threads " + THREADS + " --input " +  args.input + " --output " +  args.output + " --refdir " +  args.refdir " + --refindex " + args.refindex + " --read1 " + args.read1 + " --read2 " + args.read2
    cwd = os.getcwd()
    #print(cwd)
    if args["sso"]:
        cmd = "mpiexec -bootstrap ssh -n " + str(N) + " -ppn " + str(PPN) + \
        " --hostfile hostfile  " + \
        " python -u fq2sortedbam.py "

    else:
            
        cmd = "mpiexec -bootstrap ssh -n " + str(N) + " -ppn " + str(PPN) + \
            " -bind-to " + BINDING + \
            " -map-by " + BINDING + \
            " --hostfile hostfile  " + \
            " python -u fq2sortedbam.py "

    #print(f"N: {N}, PPN: {PPN}")
    jstring = json.dumps(args)
    try:
        #subprocess.run([cmd, 'fq2sortedbam.py', jstring, check=True, capture_output=True, text=True)
        subprocess.run([f"{cmd} '{jstring}'"], shell=True, check=True, capture_output=False, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        print(f"Error output: {e.stderr}")

    sys.exit(0)


    #cmd = "export LD_PRELOAD=" + cwd + "/libmimalloc.so.2.0:$LD_PRELOAD" + \
    cmd = "mpiexec -bootstrap ssh -n " + str(N) + " -ppn " + str(PPN) + \
        " -bind-to " + BINDING + \
        " -map-by " + BINDING + \
        " --hostfile hostfile  " + \
        " python -u fq2sortedbam.py --cpus " + str(CPUS) + \
        " --threads " + str(THREADS) + \
        " --input " +  args["input"] + \
        " --output " + args["output"] + \
        " --tempdir " +  args["tempdir"] + \
        " --refdir " +  args["refdir"] + \
        " --refindex " + args["refindex"] + \
        " --read1 " + args["read1"] + \
        " --read2 " + args["read2"] + \
        " --read3 " + args["read3"] + \
        " --readi1 " + args["readi1"] + \
        " --params " + "\"" + args["params"] + "\"" + \
        " --prefix " + args["prefix"] + \
        " --suffix " + args["suffix"] + \
        " --whitelist " + args["whitelist"] + \
        " --read_structure " + args["read_structure"] + \
        " --barcode_orientation " + args["barcode_orientation"] + \
        " --sample_id " + str(args["sample_id"]) + \
        " --output_format " + args["output_format"] + \
        " --bam_size " + str(args["bam_size"]) + \
        " --mode " + args["mode"] + \
        " --read_type " + args["read_type"] + \
        " --outfile " + args["outfile"]
        ##" --container_tool " + args["container_tool #+ \

    if args["rindex"]:
        cmd += " --rindex "
    if args["dindex"]:
        cmd +=" --dindex " 
    if not args["not_keep_unmapped"]:
        cmd += " --keep_unmapped "

    cmd += " 2>&1 | tee " + args["output"] + "/log_fq2sortedbam.txt"        
    
    print("cmd: ", cmd)
    a = run(cmd, capture_output=True, shell=True)
    assert a.returncode == 0, "fq2bam failed."
