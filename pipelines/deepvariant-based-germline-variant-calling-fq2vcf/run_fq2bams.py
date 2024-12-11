from subprocess import Popen, PIPE, run
import json, os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
#from mpi4py import MPI

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

    if sso == 'sso':
        nsocks = 1

    #th = 16  ## max cores per rank
    #num_nodes = 1
        
    num_physical_cores_all_nodes = num_nodes * nsocks * ncores
    num_physical_cores_per_node = nsocks * ncores
    num_physical_cores_per_rank = nsocks * ncores
    
    while num_physical_cores_per_rank > int(th):
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
    #parser.add_argument('--tempdir',default="/output",help="Intermediate data directory")
    #parser.add_argument('--refdir',default="/refdir",help="Reference genome directory")
    parser.add_argument('--output',default="/output", help="Output data directory")
    parser.add_argument("-i", "--refindex", default="None", help="name of refindex file")
    #parser.add_argument("-r1", "--read1", default="None", help="name of read1")
    #parser.add_argument("-r2", "--read2", default="None", help="name of read2")
    parser.add_argument('-in', '--rindex',action='store_true',help="It will index reference genome for bwa-mem2. If it is already done offline then don't use this flag.")
    parser.add_argument('-dindex',action='store_true',help="It will create .fai index. If it is done offline then disable this.")
    parser.add_argument('--container_tool',default="docker",help="Container tool used in pipeline : Docker/Podman")
    parser.add_argument('-pr', '--profile',action='store_true',help="Use profiling")
    parser.add_argument('--keep_unmapped',action='store_true',help="To not keep the unmapped entries at the end of sam file.")
    parser.add_argument('--keep_intermediate_sam',action='store_true',help="Keep intermediate sam files.")
    parser.add_argument('--params', type=str, default='@RG\\tID:RG1\\tSM:RGSN1', help="parameter string to bwa-mem2 barring threads paramter")
    #parser.add_argument("-p", "--outfile", help="prefix for read files")
    parser.add_argument("--sso", action='store_true', help="prefix for read files")
    parser.add_argument('--buildindexonly',action='store_true',help="It will create bwa and .fai index only. If it is done offline then disable this.")
    parser.add_argument("--th", default=20, help="#min cores per rank")
    parser.add_argument("-N", default=-1, help="#ranks")
    parser.add_argument("-PPN", default=-1, help="ppn")
    parser.add_argument("--cpus", default=-1, help="CPUS")
    
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
    #N, PPN, CPUS, THREADS, mask = HWConfigure(args["sso"], num_nodes)
    N, PPN, CPUS, THREADS, mask, numa_per_sock = HWConfigure(args["sso"], num_nodes, args['th'])    
    if args["N"] != -1: N = args["N"]
    if args["PPN"] != -1: PPN = args["PPN"]
    if args["cpus"] != -1: CPUS = args["cpus"]
    #if args["threads"] != -1: THREADS = args["cpus"]
    
    print("[Info] Running {} processes per compute node, each with {} threads".format(N, THREADS))
    args['cpus'], args['threads'] = str(CPUS), str(THREADS)
    
    #print("[Info] Running {} processes per compute node, each with {} threads".format(N, THREADS))
    cmd="hostname > hostfile"
    a = run(cmd, capture_output=True, shell=True)
    
    BINDING="socket"
    cmd="mkdir -p logs"
    a = run(cmd, capture_output=True, shell=True)

    lpath="/app/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/libmimalloc.so.2.0"
        
    if args["sso"]:
        print(f'Running on single socket w/ {numa_per_sock} numas per socket')
        cmd = "export LD_PRELOAD=" + lpath + "; numactl -N " + "0-" + str(numa_per_sock-1) + " mpiexec -bootstrap ssh -n " + str(N) + " -ppn " + str(PPN) + \
        " --hostfile hostfile  " + \
        " python -u fq2sortedbam.py "

    else:            
        cmd = "export LD_PRELOAD=" + lpath + "; mpiexec -bootstrap ssh -n " + str(N) + " -ppn " + str(PPN) + \
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
    
    #cmd = "export I_MPI_PIN_DOMAIN==mask" + "; mpiexec -bootstrap ssh -n " + N + "-ppn " + PPN + " -bind-to " + BINDING + "-map-by " + BINDING + " --hostfile hostfile  python -u fq2bams.py --cpus" + CPUS + " --threads " + THREADS + " --input " +  args.input + " --output " +  args.output + " --refdir " +  args.refdir " + --refindex " + args.refindex + " --read1 " + args.read1 + " --read2 " + args.read2
    cwd = os.getcwd()
    cmd = "export LD_PRELOAD=" + cwd + "/libmimalloc.so.2.0:$LD_PRELOAD" + \
        "; mpiexec -bootstrap ssh -n " + str(N) + " -ppn " + str(PPN) + \
        " -bind-to " + BINDING + \
        " -map-by " + BINDING + \
        " --hostfile hostfile  " + \
        " python -u fq2bams.py --cpus " + str(CPUS) + \
        " --threads " + str(THREADS) + \
        " --input " +  args["input"] + \
        " --output " + args["output"] + \
        " --refdir " +  args["refdir"] + \
        " --tempdir " +  args["output"] + \
        " --refindex " + args["refindex"] + \
        " --read1 " + args["read1"] + \
        " --read2 " + args["read2"] + \
        " --params " + "\"" + args["params"] + "\"" + \
        " --container_tool " + args["container_tool"] #+ \
        #" --buildindexonly " + args["buildindexonly"] #+ \
    if args["rindex"]:
        cmd += " --rindex "
    if args["dindex"]:
        cmd +=" --dindex "
    if not args["not_keep_unmapped"]:
        cmd += " --keep_unmapped "
    if args["keep_intermediate_sam"]:
        cmd +=" --keep_intermediate_sam "

    cmd += " 2>&1 | tee " + args["output"] + "/log_fq2bams.txt"        
    
    print("cmd: ", cmd)
    a = run(cmd, capture_output=True, shell=True)
    assert a.returncode == 0, "fq2bams failed."
