from subprocess import Popen, PIPE, run
import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
#from mpi4py import MPI

def HWconfigure(sso, num_nodes):
    #print("args: ", sys.argv)
    #print(len(sys.argv))
    assert len(sys.argv) == 3, "<exec> <sso> <num_nodes>"
    run('lscpu > lscpu.txt', capture_output=True, shell=True)
    #inf = sys.argv[1]
    sso = sys.argv[1]
    num_nodes = int(sys.argv[2])

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

    th = 16  ## max cores per rank
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
    SHARDS = int(ncores * nsocks / PPN)
    THREADS = CPUS
    print(f"N={int(N)}")
    print(f"PPN={int(PPN)}")
    print(f"CPUS={int(CPUS)}")
    print(f"SHARDS={int(SHARDS)}")
    print(f"THREADS={int(THREADS)}")
    
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
    print("I_MPI_PIN_DOMAIN={}".format(mask))

    return N, PPN, CPUS, THREADS, SHARDS, mask



if __name__ == '__main__':
    ## rgs parser
    parser=ArgumentParser()
    parser.add_argument('--input', default="/input", help="Input data directory")
    parser.add_argument('--tempdir',default="/tempdir",help="Intermediate data directory")
    parser.add_argument('--refdir',default="/refdir",help="Reference genome directory")
    parser.add_argument('--output',default="/output", help="Output data directory")
    parser.add_argument("-i", "--refindex", help="name of refindex file")
    #parser.add_argument("-r1", "--read1", help="name of read1")
    #parser.add_argument("-r2", "--read2", help="name of read2")
    #parser.add_argument('-in', '--rindex',action='store_true',help="It will index reference genome for bwa-mem2. If it is already done offline then don't use this flag.")
    #parser.add_argument('-dindex',action='store_true',help="It will create .fai index. If it is done offline then disable this.")
    parser.add_argument('--container_tool',default="docker",help="Container tool used in pipeline : Docker/Podman")
    parser.add_argument('-pr', '--profile',action='store_true',help="Use profiling")
    #parser.add_argument('--keep_unmapped',action='store_true',help="Keep Unmapped entries at the end of sam file.")
    #parser.add_argument('--keep_intermediate_sam',action='store_true',help="Keep intermediate sam files.")
    #parser.add_argument('--params', default='', help="parameter string to bwa-mem2 barring threads paramter")
    parser.add_argument("-p", "--outfile", help="prefix for read files")
    args = vars(parser.parse_args())

    num_nodes=1
    N, PPN, CPUS, THREADS, SHARDS, mask = HWConfigure(args.sso, num_nodes)
    cmd="hostname > hostfile"
    a = run(cmd, capture_output=True, shell=True)
    
    BINDING=socket
    cmd="mkdir -p logs"
    a = run(cmd, capture_output=True, shell=True)

    #cmd = "export I_MPI_PIN_DOMAIN==mask" + "; mpiexec -bootstrap ssh -n " + N + "-ppn " + PPN + " -bind-to " + BINDING + "-map-by " + BINDING + " --hostfile hostfile  python -u fq2bams.py --cpus" + CPUS + " --threads " + THREADS + " --input " +  args.input + " --output " +  args.output + " --refdir " +  args.refdir " + --refindex " + args.refindex + " --read1 " + args.read1 + " --read2 " + args.read2
    cwd = os.getcwd()
    cmd = "export LD_PRELOAD=" + cwd + "/libmimalloc.so.2.0:$LD_PRELOAD" + \
        "; mpiexec -bootstrap ssh -n " + N + "-ppn " + PPN + \
        " -bind-to " + BINDING + \
        " -map-by " + BINDING + \
        " --hostfile hostfile  " + \
        " python -u bams2vcf.py --cpus" + CPUS + \
        " --threads " + THREADS + \
        " --input " +  args.input + \
        " --output " + args.output + \
        " --refdir " +  args.refdir + \
        " --refindex " + args.refindex + \
        " --shards " + SHARDS + \
        " --outfile " + args.outfile + \
        " --container_tool " + args.container_tool + \
        " 2>&1 | tee " + args.output + "/log_bams2vcf.txt"
        
    
    print("cmd: ", cmd)
    a = run(cmd, capture_output=True, shell=True)
    assert a == 0, "bams2vcf failed."
    print("[Info] Pipeline executed successfully.")
