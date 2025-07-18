#!/usr/bin/env python
# encoding: utf-8

import json
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

def allexit(comm, flg):
    comm.barrier()        
    flg = comm.bcast(flg, root=0)
    if flg: os.sys.exit(1)


#def main(argv):
def main(args):
    ifile=args["refindex"]
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

    if not os.path.isfile(os.path.join(inputdir, 'aln' + binstr + '.bam')):
        print("[Info] Missing intermediate .bam files from fq2bam part of the pipeline.")
        os.sys.exit(1)
       
    if (ifile == "" or ifile == None) or not os.path.isfile(os.path.join(refdir, ifile)):
        print("[Info] Missing reference file.")
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
        os.path.join(output, 'intermediate_results_dir'+ binstr) + \
        ' --num_shards='+ str(nproc)+ \
        ' --dry_run=false --regions "' + bin_region[i*nranks+rank]+'"'
    
    print("Deepvariant commandline: ")
    print(command)
    
    a = run('echo "'+command+'" > '+os.path.join(output, "logs", 'dvlog'+binstr+'.txt'), shell=True)
    a = run(command + " 2>&1 >> " + os.path.join(output, "logs", 'dvlog'+binstr+'.txt'), shell=True)
    assert a.returncode == 0,"[Info] Deepvariant execution failed."
    comm.barrier()
    
    bins_per_rank = 1
    flg = 0
    if rank == 0:
        cmd = 'bash merge_vcf.sh '+output +' '+str(nranks)+' '+str(bins_per_rank) + ' ' + outfile + " > " + output + "/logs/mergelog.txt"
        a = run(cmd, capture_output = True, shell = True)
        #assert a.returncode == 0,"VCF merge failed"
        if a.returncode != 0:
            flg = 1
            print("[Info] VCF merge failed.")
        end5 = time.time()
        print("\nDeepVariant runtime",end5-t0)
        #print("\nTime for the whole pipeline",end5-start0)
        for i in range(nranks):
            r = "%05d"%(i)
            p = os.path.join(output, r)
            #print('path: ', p)
            os.system('rm -rf ' + p)

    if rank == nranks - 1:
        print('[Info] Cleaning up....')
        for i in range(nranks):
            r = "%05d"%(i)
            #print(r)
            #os.system('ls -lh ' + r)
            if not args['keep_input']:
                os.system('rm -rf ' + os.path.join(inputdir, "bin_region.pkl"))
                os.system('rm -rf ' + os.path.join(inputdir, "aln"+r+".bam"))
                os.system('rm -rf ' + os.path.join(inputdir, "aln"+r+".bam.bai"))

            os.system('rm -rf '+ os.path.join(output, 'intermediate_results_dir' + r))
        print('[Info] Cleaning up done.')


    allexit(comm, flg)  ## all ranks exit if failure in rank 0 above
    
if __name__ == "__main__":
    args = json.loads(sys.argv[1])
    main(args)
