from subprocess import Popen, PIPE, run
import subprocess
import time
import os
import sys
import threading
import tempfile
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pragzip as pz
from mpi4py import MPI
import bisect
import heapq as hq
import numpy as np
import yappi
from multiprocessing import Pool
from operator import itemgetter
mydict = {"chr1": 1,"chr2": 2,"chr3": 3,"chr4": 4,"chr5": 5,"chr6": 6,"chr7": 7,"chr8": 8, "chr9": 9, "chr10": 10, "chr11": 11, "chr12": 12, "chr13": 13, "chr14":14,"chr15":15,"chr16": 16, "chr17": 17, "chr18": 18, "chr19": 19,"chr20":20,"chr21":21,"chr22":22,"chrX": 23,"chrY": 24,"chrM": 25}

BINDIR="../.."


def read_input(fp_read1, fp_read2, num_ranks):
    file_path_r1 = fp_read1
    file_path_r2 = fp_read2
    lcount_r1, lcount_r2 = 0, 0
    try:
        with open(file_path_r1, 'r') as file:
            lines_r1 = file.readlines()
            lcount_r1 = len(lines_r1)
            if (lcount_r1 % num_ranks != 0):
                print("Num reads files {} is not multiple of \
                num_ranks {}".format(lcount_r1, num_ranks))

        if file_path_r2 != None:
            with open(file_path_r2, 'r') as file:
                lines_r2 = file.readlines()
                lcount_r2 = len(lines_r2)

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

    return lines_r1, lines_r2, lcount_r1, lcount_r2



def pragzip_reader_real( comm, cpus, fname, outpipe, outdir, last=False ):
    time.sleep(0.1)  # yield so bwamem can start, overlapping pragzip preindex with bwamem loading reference genome
    rank = comm.Get_rank()
    nranks = comm.Get_size()
    idxname = outdir+(fname.split('/')[-1])+".idx"
    if (rank==0 and not last) or (rank==nranks-1 and last):  # create index
        t0 = time.time()
        oldmask = os.sched_getaffinity(0)
        maxcpus = os.cpu_count()
        allcpus = {x for x in range(maxcpus)}
        os.sched_setaffinity(0,allcpus)
        f = pz.open( fname, parallelization=maxcpus//2 )
        f.export_index(idxname)
        f.close()
        os.sched_setaffinity(0,oldmask)
        t1 = time.time()
        print("Time for pragzip index: "+str(t1-t0))
    comm.Barrier()
    #f = pz.open( fname, parallelization=cpus//2 )
    f = pz.open( fname, parallelization=cpus )
    f.import_index(idxname)
    f.seek(0, 2) # go to end
    #f.seek(133430432859) # go to end
    #f.seek(133420000000) # go to end
    total_len = f.tell()
    startpos = total_len*rank//nranks  # approx starting point for each rank
    f.seek(startpos)
    while 1:
        d = f.readline()
        if d[0] == 64:  # @ symbol, begin new record
            break
        startpos = f.tell()
    #print (rank, startpos, d[0:32])
    st_vec = comm.allgather(startpos)  # get all ranks start positions
    if rank==nranks-1:
        endpos = total_len
    else:
        endpos = st_vec[rank+1]
    bufs = [0, 0]
    sems = [threading.Semaphore(2), threading.Semaphore(0)]
    threading.Thread(target=pr_t1, args=(f,startpos,endpos,bufs,sems)).start()
    threading.Thread(target=pr_t2, args=(outpipe,bufs,sems)).start()

def pragzip_reader(comm, cpus, fname, outdir, last=False):
    #outpipe = tempfile.mktemp()
    outpipe = tempfile.NamedTemporaryFile()
    temp=outpipe.name
    outpipe.close()
    #os.unlink(outpipe.name)
    outpipe=temp
    os.mkfifo(outpipe)
    comm2 = comm.Clone()
    threading.Thread(target=pragzip_reader_real, args=(comm2, cpus, fname, outpipe, outdir, last)).start()
    return outpipe



def main(argv):
    parser=ArgumentParser()
    parser.add_argument('--input',help="Input data directory")
    parser.add_argument('--temp',default="",help="Intermediate data directory")
    parser.add_argument('--refdir',default="",help="Reference genome directory")
    parser.add_argument('--output',help="Output data directory")
    parser.add_argument('--mode', str='', help="flatmode/fqprocess/pragzip. flatmode is just bwa w/o sort.")
    parser.add_argument("-i", "--index", help="name of index file")
    parser.add_argument("-p", "--rprefix", help="prefix for read files")
    parser.add_argument("-r", "--reads", nargs='+',help="name of reads file seperated by space")
    parser.add_argument("-c", "--cpus",default=1,help="Number of cpus. default=1")
    parser.add_argument("-t", "--threads",default=1,help="Number of threads used in samtool operations. default=1")
    parser.add_argument('-in', '--istart',action='store_true',help="It will index reference genome for bwa-mem2. If it is already done offline then don't use this flag.")
    #parser.add_argument('-sindex',action='store_true',help="It will create .fai index. If it is done offline then disable this.")
    #parser.add_argument('--container_tool',default="docker",help="Container tool used in pipeline : Docker/Podman")
    #parser.add_argument('--shards',default=1,help="Number of shards for deepvariant")
    parser.add_argument('-pr', '--profile',action='store_true',help="Use profiling")
    parser.add_argument('--keep_unmapped',action='store_true',help="Keep Unmapped entries at the end of sam file.")
    args = vars(parser.parse_args())
    ifile=args["index"]
    rfile1=args["reads"][0]
    rfile2=args["reads"][1]
    rprefix=args["rprefix"]
    cpus=args["cpus"]
    threads=args["threads"]    ## read prefix for R1, I1, R2 files
    istart=args["istart"]
    #sindex=args["sindex"]
    #nproc=args["shards"]
    folder=args["input"]
    output=args["output"]
    tempdir=args["temp"]
    if tempdir=="": tempdir=output
    refdir=args["refdir"]
    if refdir=="": refdir=folder
    prof=args["profile"]
    global keep
    keep=args["keep_unmapped"]
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()

    #global bins_per_rank
    #bins_per_rank = max(4,bins_per_rank//nranks)
    global ncpus
    ncpus = int(cpus)

    start0 = time.time()

     if rank==0:
         yappi.set_clock_type("wall")
         if prof: yappi.start()
         file_size = os.path.getsize(folder+rfile1)
         print("\nSize of FASTQ file:",file_size)

         if istart==True :
             print("Indexing Starts")
             begin = time.time()
             a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 index '+refdir+ifile,capture_output=True,shell=True)
             end=time.time()
             file_size = os.path.getsize(folder+rfile1)
             print("\nIndex time:",end-begin)
             aprint("\nSize of FASTQ file:",file_size)

    if mode == 'flatmode':
        if rank == 0:
            print("dist bwa-mem2 starts in flatmode")

        # Execute bwamem2 -- may include sort, merge depending on mode
        begin0 = time.time()
        fn1 = pragzip_reader( comm, int(cpus), folder+rfile1, output )
        fn2 = pragzip_reader( comm, int(cpus), folder+rfile2, output, last=True )
        fn3 = os.path.join(output, aln + str("%05d"rank) + ".sam")
        begin1 = time.time()
        #a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3,capture_output=True, shell=True)
        #bwastr = '{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'/bwalog' + str(rank) + '.txt'
        #print(bwastr)
        a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)
        assert a.returncode == 0
        end1b=time.time()
        #thr.join()
        comm.barrier()
        end1=time.time()

        if rank==0:
            print("\nFASTQ to SAM time (flatmode):",end1-begin1)
            print("   (includes wait time:",end1-end1b,")")

            print("\nsam to sort-bam starts")
            begin2=time.time()

        return

    elif mode == 'fqprocess':
        ## fastq preprocess and split
        if rank == 0:
            f'{BINDIR}/applications/warp-tools/bin/fastqprocess --verbose --bam-size '+  str(bam_size) + ' --barcode-length 16  --umi-length 10 --sample-id L8TX  --I1 ' \
                + folder/rprefix + '_I1.fastq.gz ' +  ' --R1 ' + folder/rprefix + '_R1.fastq.gz ' +   ' --R2 ' + folder/rprefix + '_R2.fastq.gz'   ## fix this based on inputs from Aseel
            ## Create #output files equal to the number of ranks
            ## create output file names as rprefix + str(rank) + '_1/2.fastq.gz'
        # Preindex refernce genome if requested
        if rank==0:
            print("bwa-mem2 starts")

        # Execute bwamem2 -- may include sort, merge depending on mode
        begin0 = time.time()
        #fn1 = pragzip_reader( comm, int(cpus), folder+rfile1, output )
        #fn2 = pragzip_reader( comm, int(cpus), folder+rfile2, output, last=True )
        fn1 = os.path.join(folder, rprefix + str(rank) + "_1.fastq")
        fn2 = os.path.join(folder, rprefix + str(rank) + "_2.fastq")
        assert os.path.isfile(fn1) == True
        assert os.path.isfile(fn2) == True
        fn3, thr = sam_writer( comm, output+'/aln' )
        begin1 = time.time()
        #a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3,capture_output=True, shell=True)
        #bwastr = '{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'/bwalog' + str(rank) + '.txt'
        #print(bwastr)
        a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)
        end1b=time.time()
        assert a.returncode == 0

        thr.join()
        comm.barrier()
        end1=time.time()

        if rank==0:
            print("\nFASTQ to SAM time (fqsplit):",end1-begin1)
            print("   (includes wait time:",end1-end1b,")")

            print("\nsam to sort-bam starts")
            begin2=time.time()

    else:
        # Preindex refernce genome if requested

        if rank == 0:
            print("bwa-mem2 starts")

        # Execute bwamem2 -- may include sort, merge depending on mode
        begin0 = time.time()
        fn1 = pragzip_reader( comm, int(cpus), folder+rfile1, output )
        fn2 = pragzip_reader( comm, int(cpus), folder+rfile2, output, last=True )
        fn3, thr = sam_writer( comm, output+'/aln' )
        begin1 = time.time()
        #a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3,capture_output=True, shell=True)
        #bwastr = '{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'/bwalog' + str(rank) + '.txt'
        #print(bwastr)
        a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)
        assert a.returncode == 0
        end1b=time.time()
        thr.join()
        comm.barrier()
        end1=time.time()

        if rank==0:
            print("\nFASTQ to SAM time:",end1-begin1)
            print("   (includes wait time:",end1-end1b,")")

            print("\nsam to sort-bam starts")
            begin2=time.time()


    # Finish sort, merge, convert to bam depending on mode
    cmd=""
    for i in range(bins_per_rank):
        binstr = '%05d'%(nranks*i+rank)
        cmd+=f'{BINDIR}/applications/samtools/samtools sort --threads '+threads+' -T '+tempdir+'/aln'+binstr+'.sorted -o '+ output +'/aln'+binstr+'.bam '+ output+'/aln'+binstr+'.sam;'
        if i%20==0:
            a=run(cmd,capture_output=True,shell=True)
            cmd=""
    if not cmd=="": a=run(cmd,capture_output=True,shell=True)
    comm.barrier()

    if rank==0:
        end2=time.time()
        print("SAM to sort-BAM time:",end2-begin2)

    ## concat bams
    if rank == 0:
        bf = []
        print('Mergining the bam files, TBD')
        for b in range(bins_per_rank):
            for r in range(nranks):
                binstr = '%05d'%(nranks*b + r)
                bf.append(output+'/aln'+binstr+'.sam')

        infstr = bf[0]
        for i in range(1, len(bf)):
            infstr = infstr + " " + bf[i]
        cmd+=f'{BINDIR}/applications/samtools/samtools cat -o' + os.path.join(output,prefix) + '.sorted.bam ' + infstr




def concatenate_files(input_files, output_file):
    try:
        with open(output_file, 'wb') as output:
            for input_file in input_files:
                with open(input_file, 'rb') as file:
                    data = file.read()
                    output.write(data)
                # Optionally, you can insert a separator (e.g., newline) between files.
                output.write(b'\n')  # Add a newline between concatenated files

        print(f"Concatenated {len(input_files)} files into '{output_file}'.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    main(sys.argv[1:])
