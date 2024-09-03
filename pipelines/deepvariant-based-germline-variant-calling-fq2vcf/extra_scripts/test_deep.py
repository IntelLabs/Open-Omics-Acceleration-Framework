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
    parser.add_argument('--inputdir',default="",help="input to the bam folder")
    parser.add_argument('--refdir',default="",help="Reference genome directory")
    parser.add_argument('--output',help="Output data directory")
    parser.add_argument("-i", "--index", help="name of index file")
    parser.add_argument("-r", "--reads", nargs='+',help="name of reads file seperated by space")
    parser.add_argument("-c", "--cpus",default=1,help="Number of cpus. default=1")
    parser.add_argument("-t", "--threads",default=1,help="Number of threads used in samtool operations. default=1")
    parser.add_argument('--shards',default=1,help="Number of shards for deepvariant")
    args = vars(parser.parse_args())
    ifile=args["index"]
    rfile1=args["reads"]

    cpus=args["cpus"]
    threads=args["threads"]
    #istart=args["istart"]
    #sindex=args["sindex"]
    nproc=args["shards"]
    inputdir=args["inputdir"]
    output=args["output"]
    refdir=args["refdir"]
    tempdir=output

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()
    bin_region=None

    global ncpus
    ncpus = int(cpus)
    i=0
    binstr = "%05d"%(i*nranks+rank)
    t0=time.time()
    with open(os.path.join(output,'bin_region.pkl'), 'rb') as f:
        bin_region = pickle.load(f)
    print(bin_region)

    command='mkdir -p '+os.path.join(output,binstr)+'; /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref='+refdir+'/'+ifile+' --reads='+output+'/aln'+binstr+'.bam --output_vcf='+output+'/'+binstr+'/output.vcf.gz --intermediate_results_dir '+tempdir+'/intermediate_results_dir'+binstr+' --num_shards='+nproc+' --dry_run=false --regions "'+bin_region[i*nranks+rank]+'"'
    print(command)
    a = run( 'echo "'+command+'" > '+os.path.join(output,'log'+binstr+'.txt'), shell=True)
    a = run( command+" 2>&1 >> "+os.path.join(output,'log'+binstr+'.txt'), shell=True)
    assert a.returncode==0,"Deepvariant Failed"
    comm.barrier()
    bins_per_rank=1
    if rank==0:
        cmd= 'bash merge_vcf.sh '+output +' '+str(nranks)+' '+str(bins_per_rank)
        a= run(cmd,capture_output=True,shell=True)
        assert a.returncode==0,"VCF merge failed"
        end5=time.time()
        print("\nDeepVariant runtime",end5-t0)
        #print("\nTime for the whole pipeline",end5-start0)
if __name__ == "__main__":
    main(sys.argv[1:])
