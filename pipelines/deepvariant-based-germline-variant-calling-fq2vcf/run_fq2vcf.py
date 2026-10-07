#*************************************************************************************
#                           The MIT License
#
#   Intel OpenOmics - fq2vcf pipeline (unified multi-aligner, single-container)
#   Copyright (C) 2023  Intel Corporation.
#
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the
#   "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to
#   permit persons to whom the Software is furnished to do so, subject to
#   the following conditions:
#
#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#*****************************************************************************************/

import subprocess
from subprocess import Popen, PIPE, run
import json, os, sys, time
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

    num_physical_cores_all_nodes = num_nodes * nsocks * ncores
    num_physical_cores_per_node = nsocks * ncores
    num_physical_cores_per_rank = nsocks * ncores

    th = int(th)
    if th > num_physical_cores_per_rank:
        th = num_physical_cores_per_rank
        print("Threshold setting > #cores, re-setting threshold to ", th, "(num_physical_cores_per_rank)")

    while num_physical_cores_per_rank > th:
        num_physical_cores_per_rank /= 2

    num_physical_cores_per_rank = int(num_physical_cores_per_rank)
    assert num_physical_cores_per_rank > 8, 'cores per rank should be > 8'

    N = int(num_physical_cores_all_nodes / num_physical_cores_per_rank)
    PPN = int(num_physical_cores_per_node / num_physical_cores_per_rank)
    CPUS = int(ncores * nthreads * nsocks / PPN - 2*nthreads)
    THREADS = CPUS

    threads_per_rank = num_physical_cores_per_rank * nthreads
    bits = pow(2, num_physical_cores_per_rank) - 1
    allbits = 0
    mask="["
    for r in range(N):
        allbits = allbits | (bits << r*num_physical_cores_per_rank)
        allbits = allbits | (allbits << nsocks * ncores)
        if mask == "[":
            mask = mask + hex(allbits)
        else:
            mask = mask+","+ hex(allbits)
        allbits=0
    mask=mask + "]"

    return N, PPN, CPUS, THREADS, mask, numa_per_sock


if __name__ == '__main__':
    parser=ArgumentParser()
    parser.add_argument('--ref', default="", help="Reference genome path. For bwa-mem2/bwa-meth this pipeline expects the index here (use --rindex to build it). For mm2-fast/STAR, pre-build the index/genomeDir separately.")
    parser.add_argument('--reads', nargs='+', help="Input reads, expects both the reads at the same location.")
    parser.add_argument('--output',default="/output/out.vcf", help="Output vcf (or bam, for STAR/bwa-meth) prefix location.")
    parser.add_argument('--simd',default="avx", help="Defaults to avx512 mode, use 'sse' for bwa-mem2 sse mode.")
    parser.add_argument("--refindex", default="None", help="name of refindex file")
    parser.add_argument('--aligner', default="bwa-mem2", choices=["bwa-mem2", "mm2-fast", "bwa-meth", "STAR"],
                         help="Aligner to use. bwa-mem2/mm2-fast feed into DeepVariant for variant calling (WGS/PACBIO model respectively); bwa-meth/STAR output a sorted BAM only.")
    parser.add_argument('--model_type', default="", help="Override DeepVariant --model_type (only used when --aligner is bwa-mem2 or mm2-fast). Defaults to WGS for bwa-mem2, PACBIO for mm2-fast.")
    parser.add_argument('--rindex',action='store_true',help="Builds the reference index for bwa-mem2/bwa-meth. Not supported for mm2-fast/STAR -- pre-build their indexes separately.")
    parser.add_argument('-dindex',action='store_true',help="It will create .fai index. If it is done offline then disable this.")
    parser.add_argument('--profile',action='store_true',help="Use profiling")
    parser.add_argument('--not_keep_unmapped',action='store_true',help="It rejects unmapped reads at the end of sorted bam file, else it accepts the unmapped reads.")
    parser.add_argument('--keep_intermediate_sam',action='store_true',help="It keeps intermediate SAM files generated out of the alignment tool for each rank. SAM file naming: aln{rank:04d}.sam")
    parser.add_argument('--keep_input', action='store_true', help="Keep intermediate per-rank BAM files after DeepVariant calling. By default they are removed.")
    parser.add_argument('--params', type=str, default='-R "@RG\\tID:RG1\\tSM:RGSN1"', help="Enables supplying various parameters to the chosen aligner (barring threads parameter). e.g. --params '-R \"@RG\\tID:RG1\\tSM:RGSN1\"\' for read grouping.")
    parser.add_argument("--sso", action='store_true', help="Uses only single socket for execution.")
    parser.add_argument('--buildindexonly',action='store_true',help="It will create bwa-mem2 and .fai index only. If it is done offline then disable this.")
    parser.add_argument("--th", default=20, help="Threshold for minimum cores allocation to each rank")
    parser.add_argument("-N", default=-1, help="Enables manual setting of #ranks. While using this setting please set PPN, cpus options accordingly.")
    parser.add_argument("-PPN", default=-1, help="Enables manual setting of ppn. While using this setting please set N, cpus options accordingly.")
    parser.add_argument("--cpus", default=-1, help="Enables manual setting of cpus option. While using this setting please set N, PPN options accordingly.")

    args = vars(parser.parse_args())

    assert len(args["reads"]) >= 1

    args["input"] = os.path.dirname(args["reads"][0])

    args["read1"] = os.path.basename(args["reads"][0])
    if len(args["reads"]) == 2: args["read2"] = os.path.basename(args["reads"][1])
    else: args["read2"] = ""

    args["refdir"] = os.path.dirname(args["ref"])
    args["refindex"] = os.path.basename(args["ref"])
    args["output"] = os.path.dirname(args["output"])
    args["tempdir"] = args["output"]
    args["outfile"] = os.path.basename(args["output"])

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

    lpath="/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/libmimalloc.so.2.0"

    if args["sso"]:
        print(f'Running on single socket w/ {numa_per_sock} numas per socket')
        cmd = "export LD_PRELOAD=" + lpath + "; numactl -N " + "0-" + str(numa_per_sock-1) + " mpiexec -bootstrap ssh -n " + str(N) + " -ppn " + str(PPN) + \
        " --hostfile hostfile  " + \
        " python -u fq2vcf.py "
    else:
        cmd = "export LD_PRELOAD=" + lpath + "; mpiexec -bootstrap ssh -n " + str(N) + " -ppn " + str(PPN) + \
            " -bind-to " + BINDING + \
            " -map-by " + BINDING + \
            " --hostfile hostfile  " + \
            " python -u fq2vcf.py "

    jstring = json.dumps(args)
    tic = time.time()
    try:
        subprocess.run([f"{cmd} '{jstring}'"], shell=True, check=True, capture_output=False, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        print(f"Error output: {e.stderr}")
    toc = time.time()

    print('[Info] fq2vcf runtime: {:.2f}'.format(toc - tic), " sec")
    sys.exit(0)
