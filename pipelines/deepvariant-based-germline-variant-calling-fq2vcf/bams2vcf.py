#!/usr/bin/env python
# encoding: utf-8

from subprocess import Popen, PIPE, run
import subprocess
import time
import os
import sys
import threading
import tempfile
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from mpi4py import MPI
import bisect
import heapq as hq
import numpy as np
from multiprocessing import Pool
from operator import itemgetter
import pickle
BINDIR="../.."

def main(argv):
    parser=ArgumentParser()
    parser.add_argument('--input', default="/input", help="Input data directory")
    parser.add_argument('--tempdir',default="/tempdir",help="Intermediate data directory")
    parser.add_argument('--refdir',default="/refdir",help="Reference genome directory")
    parser.add_argument('--output',default="/output", help="Output data directory")
    parser.add_argument("-i", "--refindex", help="name of refindex file")
    #parser.add_argument("-r", "--reads", nargs='+',help="name of reads file seperated by space")
    #parser.add_argument("-c", "--cpus",default=1,help="Number of cpus. default=1")
    #parser.add_argument("-t", "--threads",default=1,help="Number of threads used in samtool operations. default=1")
    #parser.add_argument('--shards',default=1,help="Number of shards for deepvariant")
    parser.add_argument("-p", "--outfile", help="prefix for the output vcf file")    
    args = vars(parser.parse_args())
    ifile=args["refindex"]
    #rfile1=args["reads"]

    cpus=args["cpus"]
    threads=args["threads"]
    nproc=args["shards"]
    inputdir=args["input"] + "/"
    output=args["output"] + "/"
    refdir=args["refdir"] + "/"
    #tempdir=output
    tempdir = args["tempdir"]
    if tempdir == "": tempdir = output
    else: tempdir = tempdir + "/"

    outfile = args["outfile"]

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()
    bin_region=None

    global ncpus
    ncpus = int(cpus)
    i = 0
    binstr = "%05d"%(i*nranks+rank)
    t0 = time.time()
    if not os.path.isfile(os.path.join(output, 'bin_region.pkl')):
        print("[Info] Missing intermediate .pkl files from fq2bam part of the pipeline.")
        os.sys.exit(1)

    if not os.path.isfile(os.path.join(inputdir, 'aln' + binstr + '.bam'))
        print("[Info] Missing intermediate .bam files from fq2bam part of the pipeline.")
        os.sys.exit(1)
        
    with open(os.path.join(inputdir, 'bin_region.pkl'), 'rb') as f:
        bin_region = pickle.load(f)
    print(bin_region)

    command='mkdir -p '+os.path.join(output,binstr)+ \
        '; /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=' + \
        os.path.join(refdir, ifile) + \
        ' --reads='+ os.path.join(inputdir, 'aln' + binstr + '.bam') + ' '  + \
        ' --output_vcf=' + os.path.join(output, binstr, 'output.vcf.gz ') + ' ' + \
        ' --intermediate_results_dir '+ \
        os.path.join(tempdir, '/intermediate_results_dir'+ binstr) + \
        ' --num_shards='+nproc+ \
        ' --dry_run=false --regions "' + bin_region[i*nranks+rank]+'"'
    
    print("Deepvariant commandline: ")
    print(command)
    
    a = run('echo "'+command+'" > '+os.path.join(output, 'dvlog'+binstr+'.txt'), shell=True)
    a = run(command + " 2>&1 >> " + os.path.join(output, 'dvlog'+binstr+'.txt'), shell=True)
    assert a.returncode == 0,"[Info] Deepvariant execution failed."
    comm.barrier()
    
    bins_per_rank = 1
    flg = 0
    if rank == 0:
        cmd = 'bash merge_vcf.sh '+output +' '+str(nranks)+' '+str(bins_per_rank)
        a = run(cmd, capture_output = True, shell = True)
        #assert a.returncode == 0,"VCF merge failed"
        if a.returncode != 0:
            flg = 1
            print("[Info] VCF merge failed.")
        end5 = time.time()
        print("\nDeepVariant runtime",end5-t0)
        #print("\nTime for the whole pipeline",end5-start0)

    allexit(comm, flg)  ## all ranks exit if failure in rank 0 above

    
if __name__ == "__main__":
    main(sys.argv[1:])
