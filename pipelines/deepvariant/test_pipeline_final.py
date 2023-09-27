#!/usr/bin/python
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
    time.sleep(0.1)  # yield so bwamem ca start, overlapping pragzip preindex with bwamem loading reference genome
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
    parser.add_argument('--output',help="Output data directory")
    parser.add_argument("-i", "--index", help="name of index file")
    parser.add_argument("-r", "--reads", nargs='+',help="name of reads file seperated by space")
    parser.add_argument("-c", "--cpus",default=1,help="Number of cpus. default=1")
    parser.add_argument("-t", "--threads",default=1,help="Number of threads used in samtool operations. default=1")
    parser.add_argument('-in', '--istart',action='store_true',help="It will index reference genome for bwa-mem2. If it is already done offline then don't use this flag.") 
    parser.add_argument('-sindex',action='store_true',help="It will create .fai index. If it is done offline then disable this.")
    parser.add_argument('--container_tool',default="docker",help="Container tool used in pipeline : Docker/Podman")
    parser.add_argument('--shards',default=1,help="Number of shards for deepvariant")
    parser.add_argument('-pr', '--profile',action='store_true',help="Use profiling") 
    parser.add_argument('--keep_unmapped',action='store_true',help="Keep Unmapped entries at the end of sam file.")
    args = vars(parser.parse_args())
    ifile=args["index"]
    rfile1=args["reads"][0]
    rfile2=args["reads"][1]
    cpus=args["cpus"]
    threads=args["threads"]
    istart=args["istart"]
    sindex=args["sindex"]
    nproc=args["shards"]
    folder=args["input"]
    output=args["output"]
    tempdir=args["temp"]
    if tempdir=="": tempdir=output
    refdir=args["refdir"]
    if refdir=="": refdir=folder
    prof=args["profile"]
    global keep
    keep=args["keep_unmapped"]
    container_tool=args["container_tool"]
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()

    #global bins_per_rank
    #bins_per_rank = max(4,bins_per_rank//nranks)
    global ncpus
    ncpus = int(cpus)

    start0 = time.time()

    # Preindex refernce genome if requested
    if rank==0:
        yappi.set_clock_type("wall")
        if prof: yappi.start()
        file_size = os.path.getsize(os.path.join(folder,rfile1))
        print("\nSize of FASTQ file:",file_size)
        
        if istart==True :
            print("Indexing Starts")
            begin = time.time()
            a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 index '+os.path.join(refdir,ifile),capture_output=True,shell=True)
            end=time.time()
            file_size = os.path.getsize(os.path.join(folder,rfile1))
            print("\nIndex time:",end-begin)
            print("\nSize of FASTQ file:",file_size)
            # numactl -m 0 -N 0 ./bwa-mem2/bwa-mem2 index data/tempdata/refs/hs37d5.fa.gz
        
        print("bwa-mem2 starts")


    # Execute bwamem2 -- may include sort, merge depending on mode
    begin0 = time.time()
    fn1 = pragzip_reader( comm, int(cpus), os.path.join(folder,rfile1), output )
    fn2 = pragzip_reader( comm, int(cpus), os.path.join(folder,rfile2), output, last=True )
    fn3, thr = sam_writer( comm, os.path.join(tempdir,'aln') )
    begin1 = time.time()
    a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+os.path.join(refdir,ifile)+' '+fn1+' '+fn2+' > '+fn3,capture_output=True, shell=True)
    end1b=time.time()
    #print("rank: ",rank," bwa_time: ",end1b-begin1)
    thr.join()
    comm.barrier()
    end1=time.time()
    #print("rank: ",rank," wait_time: ",end1-end1b)

    if rank==0:
        #file_size = os.path.getsize(folder+'aln.sam')
        #if mode>0:
        #    print("\nPreindex source files time:",begin1-begin0)
        print("\nFASTQ to SAM time:",end1-begin1)
        print("   (includes wait time:",end1-end1b,")")
        
        print("\nsam to sort-bam starts")
        begin2=time.time()
        

    # Finish sort, merge, convert to bam depending on mode
    cmd=""
    for i in range(bins_per_rank):
        binstr = '%05d'%(nranks*i+rank)
        cmd+=f'{BINDIR}/applications/samtools/samtools sort --threads '+threads+' -T '+os.path.join(tempdir,'aln'+binstr+'.sorted') + ' -o '+ os.path.join(tempdir,'aln'+binstr+'.bam')+' '+ os.path.join(tempdir,'aln'+binstr+'.sam')
        if i%20==0:
            a=run(cmd,capture_output=True,shell=True)
            cmd=""
    if not cmd=="": a=run(cmd,capture_output=True,shell=True)
    comm.barrier()
    
    if rank==0:
        end2=time.time()
        print("SAM to sort-BAM time:",end2-begin2)
        

    # Generate index file(s)
    if rank==0:
        begin3=time.time()
        print("\nIndexing of ref and read starts")
        if sindex==True :
            a=run(f'{BINDIR}/applications/samtools/samtools faidx '+os.path.join(refdir,ifile),capture_output=True,shell=True)
            end=time.time()
            print("\nReference to .fai index creation time",end-begin3)
    
    if rank==0:
        begin3=time.time()
    cmd =""
    for i in range(bins_per_rank):
        fname = "aln%05d.bam"%(i*nranks+rank)
        #a=run(f'{BINDIR}/applications/samtools/samtools index -M -@ '+threads+' '+folder+fname,capture_output=True,shell=True)
        cmd+=f'{BINDIR}/applications/samtools/samtools index -M -@ '+threads+' '+os.path.join(tempdir,fname)+';'
        if i%20==0:
            a=run(cmd,capture_output=True,shell=True)
            cmd=""
    if not cmd=="": a=run(cmd,capture_output=True,shell=True)
    comm.barrier()

    if rank==0:
        end3=time.time()
        print("\nBAM to .bai index creation time",end3-begin3)
        print("\nStarting Deepvariant execution...")
        begin5=time.time()  

    
    for i in range(bins_per_rank):
        binstr = "%05d"%(i*nranks+rank)
        command='mkdir -p '+os.path.join(output,binstr)+'; '+container_tool +' run -v '+folder+':"/input" -v '+refdir+':"/refdir" -v '+output+'/'+binstr+':"/output" -v '+tempdir+':"/tempdir" deepvariant:latest /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/refdir/'+ifile+' --reads=/tempdir/aln'+binstr+'.bam --output_vcf=/output/output.vcf.gz --intermediate_results_dir /tempdir/intermediate_results_dir'+binstr+' --num_shards='+nproc+' --dry_run=false --regions "'+bin_region[i*nranks+rank]+'" --pcl_opt'
        a = run( 'echo "'+command+'" > '+output+'log'+binstr+'.txt', shell=True)
        a = run( command+" 2>&1 >> "+output+"log"+binstr+".txt", shell=True)
    comm.barrier()
    if rank==0:
        cmd= 'bash merge_vcf.sh '+output +' '+str(nranks)+' '+str(bins_per_rank)
        a= run(cmd,capture_output=True,shell=True)
        end5=time.time()
        print("\nDeepVariant runtime",end5-begin5)

        print("\nTime for the whole pipeline",end5-start0)

        if prof:
            yappi.get_func_stats().print_all(columns={0:("ncall",5),1:("tsub",8),2:("ttot",8),3:("tavg",8),4:("name",90)})  
            yappi.get_thread_stats().print_all() 
    
if __name__ == "__main__":
    main(sys.argv[1:])

