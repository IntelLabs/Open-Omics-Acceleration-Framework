import json
import subprocess
from subprocess import Popen, PIPE, run
import os, sys, time
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


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
    if th > num_physical_cores_per_rank:
        th = num_physical_cores_per_rank
        print("Threshold setting > #cores, re-setting threshold to ", th, "(num_physical_cores_per_rank)")
        
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

    return N, PPN, CPUS, THREADS, mask, numa_per_sock



if __name__ == '__main__':
    ## rgs parser
    parser=ArgumentParser()
    parser.add_argument('--ref', default="", help="Reference genome path. For BWA this pipeline expects the index here. If the index is not present then it can be generated using --rindex option.")
    parser.add_argument('--reads', nargs='+', help="Input reads, expects both the reads at the same location.")
    #parser.add_argument('--input', default="/input", help="Input data directory")
    parser.add_argument('--tempdir',default="",help="Dir for intermediate data.")
    #parser.add_argument('--refdir',default="/refdir",help="Reference genome directory")
    parser.add_argument('--output',default="/output/out.bam", help="prefix location for the output file(s) name.")
    parser.add_argument('--simd',default="avx", help="Defaults to avx512 mode, use 'sse' for bwa sse mode.")
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
                        help="Exeuction mode: flatmode/sortedbam/fqprocessonly/multifq. flatmode: just bwa w/o sorting, creates sam files equal to the number of ranks created.\
                        sortedbam: bwa + samsort steps creating single bam file as output. Ignore other (custom) modes.")
    parser.add_argument('--read_type',default="short",
                        help="(short/long/rna/meth): short - bwa-mem2 alignment for short reads; long - mm2-fast alignment for long reads.\
                        rna - STAR alignment of rnaseq reads; meth - bwa-meth (mem2) alignment of short meth reads.  Defaults to short.")
    parser.add_argument('--rindex',action='store_true',help="It enables BWA-MEM2 index generation. Use this option if the index is not present.")
    parser.add_argument('--dindex',action='store_true',help="It creates reference genome fai index. Use this option if the reference fai is not present.")
    #parser.add_argument('--container_tool',default="docker",help="Container tool used in pipeline : Docker/Podman")
    parser.add_argument('--profile', action='store_true',help="Use profiling")
    parser.add_argument('--not_keep_unmapped',action='store_true',help="It rejects unmapped reads at the end of sorted bam file, else it accepts the unmapped reads.")
    parser.add_argument('--keep_sam',action='store_true',help="It keeps intermediate SAM files generated out of the alignment tool for each rank. SAM file naming: aln{rank:04d}.sam")    
    parser.add_argument('--params', type=str, default='-R "@RG\\tID:RG1\\tSM:RGSN1"', help="Enables supplying various parameters to bwa-mem2 (barring threads (-t) paramter). e.g. --params '-R \"@RG\\tID:RG1\\tSM:RGSN1\"\' for read grouping.")
    #parser.add_argument("-p", "--outfile", default="final", help="prefix for read files")
    parser.add_argument("--sso", action='store_true', help="Executes the pipeline on single socket only. By default, it uses all the sockets.")
    parser.add_argument("--th", default=20, help="Threshold for minimum cores allocation to each rank")
    parser.add_argument("-N", default=-1, help="Enables manual setting of #ranks. While using this setting please set PPN, cpus options accordingly.")
    parser.add_argument("-PPN", default=-1, help="Enables manual setting of ppn. While using this setting please set N, cpus options accordingly.")
    parser.add_argument("--cpus", default=-1, help="Enables manual setting of cpus option. While using this setting please set N, PPN options accordingly.")
    #parser.add_argument("--threads", default=-1, help="THREADS")
    args = vars(parser.parse_args())

    assert len(args["reads"]) >= 1
    
    args["input"] = os.path.dirname(args["reads"][0])

    args["read1"] = os.path.basename(args["reads"][0])
    if len(args["reads"]) == 2: args["read2"] = os.path.basename(args["reads"][1])
    else: args["read2"] = ""

    args["refdir"] = os.path.dirname(args["ref"])
    args["refindex"] = os.path.basename(args["ref"])
    args["outfile"] = os.path.basename(args["output"])
    args["output"] = os.path.dirname(args["output"])
    args["tempdir"] = args["output"]
    
    num_nodes=1
    N, PPN, CPUS, THREADS, mask, numa_per_sock = HWConfigure(args["sso"], num_nodes, args['th'])
    if args["N"] != -1:
        N = args["N"]
        assert args["PPN"] != -1, "Please set PPN when manually setting N"
        assert args["cpus"] != -1, "Please set cpus when manually setting N"
        
    if args["PPN"] != -1:
        PPN = args["PPN"]
        assert args["N"] != -1, "Please set N when manually setting PPN"
        assert args["cpus"] != -1, "Please set cpus when manually setting PPN"
        
    if args["cpus"] != -1:
        CPUS = args["cpus"]
        THREADS = args["cpus"]
        assert args["PPN"] != -1, "Please set PPN when manually setting cpus"
        assert args["N"] != -1, "Please set N when manually setting cpus"        

    
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
    lpath="/app/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/libmimalloc.so.2.0"
    tic = time.time()
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

    toc = time.time()
    print('[Info] fq2sortedbam runtime: {:.2f}'.format(toc - tic), " sec")
