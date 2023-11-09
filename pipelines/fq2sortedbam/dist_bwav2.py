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
#from warp-tools.tools.scripts import dynamic-barcode-orientation
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

def pr_t1( f, start, end, bufs, sems ):  # thread to read from fastq.gz file
    f.seek(start)
    buf_i = 0
    k = 1
    while k>0:
        sems[0].acquire()
        k = min(1_000_000, end-f.tell())
        if k<1:
            d = None
        else:
            d = f.read(k)
        bufs[buf_i] = d
        buf_i = 1-buf_i
        sems[1].release()

def pr_t2( outpipe, bufs, sems ):   # thread to output data to pipe
    fifo = open(outpipe, 'wb')
    buf_i = 0
    while 1:
        sems[1].acquire()
        d = bufs[buf_i]
        if d is None: break
        buf_i = 1-buf_i
        sems[0].release()
        fifo.write(d)
    fifo.close()
    os.unlink(outpipe)  # delete outpipe



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


#bins_per_rank = 500  # actually -- set this to number of bins, it will get adjusted by main to something reasonable per rank
bins_per_rank = 1
pq = []   # in-memory lists to store, sort records output by bwamem2; sort key is position; (one heap per local bin)
headers = []   # header line in SAM file;  will indicate the sequence names and their lenghts
headerlen = 0
seq_start = {}  # sequence name to cumulative start position
seq_len = {}  # sequence name to length
cumlen = 0 # max position value
bins = []  # tuples of (start_pos, rank, local bin_num)
bin_region = [] # strings of form seq_id:start-end
ncpus = 8
binrounding = 1000
keep=False
headers_done = threading.Semaphore(0)

# NOTE:  Code relies on insert order of dictionary keys -- so needs python>=3.7

# Calculate bins
def calculate_bins(nranks):
    global bins
    seq_ids = np.array(list(seq_len.keys()))
    seq_lens = np.array(list(seq_len.values()))
    seq_i = 0
    bin_i = 0
    start = 0  # start of next bin in global position numbers
    for b in range(bins_per_rank):
        for r in range(nranks):
            bin_reg = ""
            end = (bin_i+1)*cumlen // (nranks*bins_per_rank)
            while seq_start[seq_ids[seq_i]]+seq_lens[seq_i] < end-binrounding:  # this seq finishes by desired end of bin
                bin_reg+=seq_ids[seq_i]+':'+str(max(0,start-seq_start[seq_ids[seq_i]]))+'-'+str(seq_lens[seq_i])+' '
                seq_i+=1
            # round off bin end
            if abs(end-(seq_start[seq_ids[seq_i]]+seq_lens[seq_i]))<=binrounding:
                end = seq_start[seq_ids[seq_i]]+seq_lens[seq_i]
            else:
                end = (end-seq_start[seq_ids[seq_i]]+binrounding//2)//binrounding*binrounding + seq_start[seq_ids[seq_i]]
            bin_reg+=seq_ids[seq_i]+':'+str(max(0,start-seq_start[seq_ids[seq_i]]))+'-'+str(min(end-seq_start[seq_ids[seq_i]],seq_lens[seq_i]))
            bins.append( (start, r, b) )
            bin_region.append( bin_reg )
            start = end
            bin_i += 1
    #print (seq_len)
    #print (bins, bin_region)
    #assert 0

# used in mode 5
def sort_thr5(fname, comm):
    rank = comm.Get_rank()
    nranks = comm.Get_size()
    ndone = 0
    t0 = time.time()
    fl = [ open(fname+("%05d.sam"%(i*nranks+rank)),"w") for i in range(bins_per_rank) ]
    headers_done.acquire()  # wait until headers are read from bwamem2
    h = "".join(headers)
    for f in fl:
        f.write(h)
    # get binned items, output to bin file
    while ndone<nranks:
        #b, vl = comm.recv()
        req = comm.irecv()
        b, vl = req.wait()
        b = b//nranks
        if len(vl)>0 and  vl[0]=="done":
            ndone+=1
            #print ("sort_thr "+str(rank)+", got done:"+str(ndone))
        else:
            fl[b].writelines(vl)
    [f.close() for f in fl]
    t1 = time.time()
    #print("sort_thr "+str(rank)+" time: "+str(t1-t0))


# used for mode 2 and above
#sent_msgs= []
def sw_thr( outpipe, comm, comm2 ):
    rank = comm.Get_rank()
    nranks = comm.Get_size()
    global headerlen
    global cumlen
    global bins
    send_bufs = [[] for _ in range(nranks*bins_per_rank)]

    def bisect_wrap(b,k):
        return bisect.bisect(b, k )
    def send_wrap(v,r,force=False, b=None):
        #global sent_msgs
        #if r==rank: return
        if b is None: b=r
        send_bufs[b].append(v)
        if force or len(send_bufs[b])>=20:
            #print("send",rank,r)
            #comm.send( (b, send_bufs[b]),r)
            req = comm.isend( (b, send_bufs[b]),r)
            req.wait()
            #if not req.Test():
            #    sent_msgs.append( (req, send_bufs[b]) )  # if send is not yet complete, save this request and data
            #    print("saved request",rank,r)
            #req.wait()
            send_bufs[b].clear()
            #send_bufs[b] = []
            #if force: print ("send forced on rank "+str(rank)+" to "+str(r))
        #sent_msgs = [ rqt for rqt in sent_msgs if not rqt[0].Test() ]
    def read_wrap(f):
        return f.readline()

    t0 = time.time()
    #seq_start['*'] = 0    # put all unmatched reads at beginning
    cumlen = 0
    global keep
    with open(outpipe,'r') as f1:
        #l = f1.readline()
        l = read_wrap(f1)
        while l and l[0]=='@': # header lines
            headers.append(l)
            headerlen += len(l)
            l = l.split()
            if l[0] == '@SQ':   # @SQ lines describe sequences in reference genome
                sn = l[1].split(':')[1]
                ln = int(l[2].split(':')[1])
                if keep:
                    seq_start[sn] = cumlen
                    seq_len[sn] = ln
                    cumlen += ln
                else:
                    if sn in mydict.keys():
                        seq_start[sn] = cumlen
                        seq_len[sn] = ln
                        cumlen += ln
                    else :
                        seq_start[sn] = -1
		#l = f1.readline()
            l = read_wrap(f1)
        # done reading headers
        if keep:
            seq_start['*'] = cumlen    # put all unmatched reads at end
        else:
            seq_start['*'] = -1

        #print(seq_start)
        calculate_bins(nranks)
        binstarts = [ b[0] for b in bins ]
        #if rank==0: print("bins", bins)
        headers_done.release()   # allow output thread to begin
        #print (binstarts)

        # read remaining lines (aligned short reads), send to appropriate bin
        i = 0
        while l:
            x = l.split()
            seq = x[2]
            offset = int(x[3])
            #print (i, seq, offset)
            if keep:
            	key = seq_start[seq] + offset
            	#b = bisect.bisect(binstarts, key) - 1
            	bn = bisect_wrap(binstarts, key) - 1
            	#if bn==0 and not (seq=='chr1' or seq=='chr2'):
            	#    print("BAD BIN", seq, seq_start[seq], offset, key, bn)
            	_, r, b = bins[bn]
            	#print (seq, offset, key, bn, r, b)
            	#comm.send( (key, b, l), r )

            	send_wrap( l, r, b=bn )
            else :
                if seq_start[seq]!= -1:
                    key = seq_start[seq] + offset
		    #b = bisect.bisect(binstarts, key) - 1
                    temp2=time.time()
                    bn = bisect_wrap(binstarts, key) - 1
                    _, r, b = bins[bn]
		    #print (seq, offset, key, bn, r, b)
		    #comm.send( (key, b, l), r )
                    send_wrap( l, r, b=bn )
            l = read_wrap(f1)
            i+=1
    # send done signal
    for r in range(nranks):
        for b in range(bins_per_rank):
            send_wrap("", r, force=True, b=b*nranks+r)
        #comm.send("done", r)
        if r!=rank: send_wrap("done", r, force=True)
    send_wrap("done", rank, force=True)  # send to self last to avoid race condition with allreduce
    t1 = time.time()
    total = comm2.allreduce(i)
    #print("sw_thr "+str(rank)+" time: "+str(t1-t0)+" "+str(i)+" "+str(total))

# used for mode 2 and above
def sam_writer( comm, fname ):
    rank = comm.Get_rank()
    nranks = comm.Get_size()
    comm1 = comm.Clone()
    comm2 = comm.Clone()
    #outpipe = tempfile.mktemp()
    outpipe = tempfile.NamedTemporaryFile()
    temp=outpipe.name
    outpipe.close()
    #os.unlink(outpipe.name)
    outpipe=temp
    os.mkfifo(outpipe)
    threading.Thread(target=sw_thr, args=(outpipe,comm1,comm2)).start()
    thr = threading.Thread(target=sort_thr5, args=(fname,comm1))
    thr.start()
    return outpipe, thr


def main(argv):
    parser=ArgumentParser()
    parser.add_argument('--input',help="Input data directory")
    parser.add_argument('--temp',default="",help="Intermediate data directory")
    parser.add_argument('--refdir',default="",help="Reference genome directory")
    parser.add_argument('--read1',default="",help="name of r1 files seperated by spaces")
    parser.add_argument('--read2',default="",help="name of r2 files seperated by spaces")
    parser.add_argument('--read3',default="",help="name of r3 files seperated by spaces")
    parser.add_argument('--whitelist',default="whitelist.txt",help="10x whitelist file")
    parser.add_argument('--read_structure',default="16C",help="read structure")
    parser.add_argument('--barcode_orientation',default="FIRST_BP_RC",help="barcode orientation")
    parser.add_argument('--output',help="Output data directory")
    parser.add_argument('--mode', default='', help="flatmode/fqprocess/pragzip. flatmode is just bwa w/o sort.")
    parser.add_argument("-i", "--index", help="name of index file")
    parser.add_argument("-p", "--rprefix", help="prefix for read files")
    parser.add_argument("-r", "--preads", nargs='+',help="name of reads file seperated by space")
    parser.add_argument("-c", "--cpus",default=1,help="Number of cpus. default=1")
    parser.add_argument("-t", "--threads",default=1,help="Number of threads used in samtool operations. default=1")
    parser.add_argument("-b", "--bam_size",default=9, help="bam file size")
    parser.add_argument('-in', '--istart',action='store_true',help="It will index reference genome for bwa-mem2. If it is already done offline then don't use this flag.")
    #parser.add_argument('-sindex',action='store_true',help="It will create .fai index. If it is done offline then disable this.")
    #parser.add_argument('--container_tool',default="docker",help="Container tool used in pipeline : Docker/Podman")
    #parser.add_argument('--shards',default=1,help="Number of shards for deepvariant")
    parser.add_argument('-pr', '--profile',action='store_true',help="Use profiling")
    parser.add_argument('--keep_unmapped',action='store_true',help="Keep Unmapped entries at the end of sam file.")
    args = vars(parser.parse_args())
    ifile=args["index"]
    if args["preads"] != None:
        ## for flatmode and pragzip mode
        rfile1=args["preads"][0]
        rfile2=args["preads"][1]
    else:
        ## for fqprocss
        read1 = args["read1"]
        read2 = args["read2"]
        read3 = args["read3"]
    whitelist=args["whitelist"]
    read_structure=args["read_structure"]
    barcode_orientation=args["barcode_orientation"]
    bam_size=args["bam_size"]
    rprefix=args["rprefix"]
    mode=args["mode"]
    cpus=args["cpus"]
    threads=args["threads"]    ## read prefix for R1, I1, R2 files
    istart=args["istart"]
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
         #file_size = os.path.getsize(folder+rfile1)
         #print("\nSize of FASTQ file:",file_size)

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
        fn3 = os.path.join(output, 'aln' + str("%05d"%rank) + ".sam")
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
            print("Starting fqprocess....", flush=True)
            print("#############################################")
            print("Whitelist: ", whitelist)
            print("read_structure: ", read_structure)
            print("barcode_orientation: ", barcode_orientation)
            print("bam_size: ", bam_size)
            print("#############################################")
            #tic = time.time()
            #rl1, rl2, rl3 = read1.split(), read2.split(), read3.split()
            #barcode_index = rl2[0]
            ##cmd='zcat ' + folder/barcode_index '| sed -n "2~4p" | shuf -n 1000 > ' + output + '/downsample.fq'
            #bash_command = '''zcat '''+  folder + "/" + barcode_index + \
            #''' | sed -n '2~4p' | shuf -n 1000 > downsample.fq'''
            #a = run(bash_command,shell=True, check=True, executable='/bin/bash')
            #assert a.returncode == 0
            #
            ## dynamic-barcode-orientation downsample.fq  whitelist best_match.txt
            #try:
            #    cmd='python3 warp-tools/tools/scripts/dynamic-barcode-orientation.py downsample.fq ' \
            #        + whitelist + ' best_match.txt'
            #    a = run(cmd, shell=True, check=True, executable='/bin/bash')
            #    assert a.returncode == 0
            #except:
            #    print("Exception in dynamic-barcode-orientation!!")
            #
            ## Read the contents of the file into a variable
            #with open('best_match.txt', 'r') as file:
            #    barcode_choice = file.read().strip()
            #
            #r1, r2, r3="--R1 ", "--R2 ", "--R3 "
            #for r in range(len(rl1)):
            #    r1+=folder + "/" + rl1[r] + ' '
            #    r2+=folder + "/" + rl2[r] + ' '
            #    r3+=folder + "/" + rl3[r] + ' '
            #
            #print(r1)
            #print(r2)
            #print(r3)
            #cmd=f'warp-tools/tools/fastqpreprocessing/bin/fastqprocess --verbose --bam-size ' \
            #    +  str(bam_size) + \
            #    ' --sample-id ' + output + ' ' + r1 + ' ' + r2 + ' ' + r3 + \
            #    ' --output-format "FASTQ" '+  '--barcode-orientation ' + barcode_choice + \
            #    ' --read-structure ' + read_structure + '  --white-list ' + whitelist
            #
            #print('fastqprocess cmd: ', cmd)
            #a=run(cmd, shell=True, check=True, executable='/bin/bash')
            #assert a.returncode == 0
            ### Create #output files equal to the number of ranks
            ### create output file names as rprefix + str(rank) + '_1/2.fastq.gz'
            #toc = time.time()
            #print("Time for fq processing: ", toc - tic)

        # Preindex reference genome if requested
        comm.barrier()
        if rank==0:
            print("\nbwa-mem2 starts..")

        # Execute bwamem2 -- may include sort, merge depending on mode
        begin0 = time.time()
        #fn1 = pragzip_reader( comm, int(cpus), folder+rfile1, output )
        #fn2 = pragzip_reader( comm, int(cpus), folder+rfile2, output, last=True )
        #fn1 = os.path.join(folder, "fastq_R1_" + str(rank) + ".fastq.gz")
        #fn2 = os.path.join(folder, "fastq_R3_" + str(rank) + ".fastq.gz")

        fn1 = "fastq_R1_" + str(rank) + ".fastq.gz"
        fn2 = "fastq_R3_" + str(rank) + ".fastq.gz"
        #fn1 = folder + "/fastq_R1_" + str(rank) + ".fastq.gz"
        #fn2 = folder + "/fastq_R3_" + str(rank) + ".fastq.gz"
        if rank == 0:
            print("Input files: ")
            print(fn1)
            print(fn2)
        #assert os.path.isfile(fn1) == True
        #assert os.path.isfile(fn2) == True
        fn3, thr = sam_writer( comm, output+'/aln' )
        begin1 = time.time()
        #a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3,capture_output=True, shell=True)
        #bwastr = '{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'/bwalog' + str(rank) + '.txt'
        #print(bwastr)
        if os.path.isfile(fn1) == True:
            a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)
            assert a.returncode == 0
        else:
            print(f"{rank} No input file for me")
        end1b=time.time()

        thr.join()
        comm.barrier()
        end1=time.time()

        if rank==0:
            print("\nFASTQ to SAM time (fqprocess):",end1-begin1)
            print("   (includes wait time:",end1 - end1b,")")

            print("\nsam to sort-bam starts")
            begin2=time.time()

    elif mode == 'pragzip':
        # Preindex refernce genome if requested

        if rank == 0:
            print("bwa-mem2 starts")

        # Execute bwamem2 -- may include sort, merge depending on mode
        begin0 = time.time()
        fn1 = pragzip_reader( comm, int(cpus), folder+rfile1, output )
        fn2 = pragzip_reader( comm, int(cpus), folder+rfile2, output, last=True )
        fn3, thr = sam_writer( comm, output+'/aln' )
        begin1 = time.time()
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

    else:
        print("Incorrect mode specified!")
        os.sys.exit(1)

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
        tic = time.time()
        bf = []
        print('Concat the bam files:')
        for b in range(bins_per_rank):
            for r in range(nranks):
                binstr = '%05d'%(nranks*b + r)
                bf.append(output+'/aln'+binstr+'.bam')

        infstr = bf[0]
        for i in range(1, len(bf)):
            infstr = infstr + " " + bf[i]
        if rprefix == None:
            rprefix = "final"

        cmd+=f'{BINDIR}/applications/samtools/samtools cat -o ' + os.path.join(output,rprefix) + '.sorted.bam ' + infstr
        #print("merge cmd: ", cmd, flush=True)
        a=run(cmd,capture_output=True,shell=True)
        assert a.returncode == 0
        print("Time for cat: ", time.time() - tic)


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
