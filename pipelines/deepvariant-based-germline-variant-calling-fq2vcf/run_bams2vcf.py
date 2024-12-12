import json, subprocess
from subprocess import Popen, PIPE, run
import os, sys, time
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
#from mpi4py import MPI

def HWConfigure(sso, num_nodes, th=20):
    run('lscpu > lscpu.txt', capture_output=True, shell=True)
    dt={}
    flg, count = 1, -1
    numa_cpu = []

    with open('lscpu.txt', 'r') as f:
        l = f.readline()
        while l:
            try:
                a,b = l.strip('\n').split(':')
                #aa, bb = a.split(' ')
                #print(a, b)
                dt[a] = b
                if a.startswith("NUMA") == True and count > 0:
                    numa_cpu.append(b.lstrip())
                
                if a.startswith("NUMA") == True and flg and b != "":
                    flg = 0
                    nnuma = int(dt['NUMA node(s)'])
                    count = int(b)

            except:
                pass
            
            l = f.readline()

    ncpus = int(dt['CPU(s)'])
    nsocks = int(dt['Socket(s)'])
    nthreads = int(dt['Thread(s) per core'])
    ncores = int(dt['Core(s) per socket'])
    nnuma = int(dt['NUMA node(s)'])
    numa_per_sock = int(count/nsocks)
    print('CPUS: ', ncpus)
    print('#sockets: ', nsocks)
    print('#threads: ', nthreads)
    print('NUMAs: ', nnuma)    

    if sso:
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
    SHARDS = int(ncores * nthreads * nsocks / PPN) 
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

    return N, PPN, CPUS, THREADS, SHARDS, mask, numa_per_sock



if __name__ == '__main__':
    ## rgs parser
    parser=ArgumentParser()
    parser.add_argument('--ref', default="", help="Input data directory")
    parser.add_argument('--input', default="/input", help="Input data directory")
    #parser.add_argument('--tempdir',default="/output",help="Intermediate data directory")
    #parser.add_argument('--refdir',default="/refdir",help="Reference genome directory")
    parser.add_argument('--output', default="/output/out.vcf", help="Output data directory")
    #parser.add_argument("-i", "--refindex", default="None", help="name of refindex file")
    #parser.add_argument('--container_tool',default="docker",help="Container tool used in pipeline : Docker/Podman")
    parser.add_argument('-pr', '--profile',action='store_true',help="Use profiling")
    #parser.add_argument("-p", "--outfile", default="finalvcf", help="prefix for the output vcf file")
    parser.add_argument("--sso", action="store_true", help="single socket execution")
    parser.add_argument("--th", default=20, help="#min cores per rank")
    parser.add_argument("-N", default=-1, help="#ranks")
    parser.add_argument("-PPN", default=-1, help="ppn")
    parser.add_argument("--cpus", default=-1, help="cpus")
    parser.add_argument("--threads", default=-1, help="threads")
    parser.add_argument("--shards", default=-1, help="shards")
    
    args = vars(parser.parse_args())

    args["refdir"] = os.path.dirname(args["ref"])
    args["refindex"] = os.path.basename(args["ref"])
    args["output"] = os.path.dirname(args["output"])
    args["tempdir"] = os.path.dirname(args["output"])
    args["outfile"] = os.path.basename(args["output"])
    print('outfile: ', args['outfile'])

    num_nodes=1
    ## N, PPN, CPUS, THREADS, SHARDS, mask = HWConfigure(args["sso"], num_nodes, args['th'])
    N, PPN, CPUS, THREADS, SHARDS, mask, numa_per_sock = HWConfigure(args["sso"], num_nodes, args['th'])    
    if args["N"] != -1: N = args["N"]
    if args["PPN"] != -1: PPN = args["PPN"]
    if args["cpus"] != -1: CPUS = args["cpus"]
    if args["shards"] != -1: SHARDS = args["cpus"]
    
    print('[Info] Running {} processes on one compute node, each with {} threads'.format(N, THREADS))
    args['cpus'], args['threads'], args['shards'] = str(CPUS), str(THREADS), str(SHARDS)

    cmd="hostname > hostfile"
    a = run(cmd, capture_output=True, shell=True)
    
    BINDING="socket"
    cmd="mkdir -p logs"
    a = run(cmd, capture_output=True, shell=True)

    lpath="/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/libmimalloc.so.2.0"
    if args["sso"]:
        print(f'Running on single socket w/ {numa_per_sock} numas per socket')
        cmd = "export LD_PRELOAD=" + lpath + "; numactl -N " + "0-" + str(numa_per_sock-1) + " mpiexec -bootstrap ssh -n " + str(N) + " -ppn " + str(PPN) + \
        " --hostfile hostfile  " + \
        " python -u bams2vcf.py "

    else:
        cmd = "export LD_PRELOAD=" + lpath + "; mpiexec -bootstrap ssh -n " + str(N) + " -ppn " + str(PPN) + \
            " -bind-to " + BINDING + \
            " -map-by " + BINDING + \
            " --hostfile hostfile  " + \
            " python -u bams2vcf.py "
    
    jstring = json.dumps(args)
    try:
        subprocess.run([f"{cmd} '{jstring}'"], shell=True, check=True, capture_output=False, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        print(f"Error output: {e.stderr}")

    print("[Info] bams2vcf pipeline executed successfully.")
