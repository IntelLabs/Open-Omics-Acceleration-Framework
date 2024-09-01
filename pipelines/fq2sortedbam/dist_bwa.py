#*************************************************************************************
#                           The MIT License
#
#   Intel OpenOmics - fq2sortedbam pipeline
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
#
#Authors:  Vasimuddin Md <vasimuddin.md@intel.com>; Babu Pillai <padmanabhan.s.pillai@intel.com>;
#*****************************************************************************************/


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
import yaml
#from warp-tools.tools.scripts import dynamic-barcode-orientation
mydict = {"chr1": 1,"chr2": 2,"chr3": 3,"chr4": 4,"chr5": 5,"chr6": 6,"chr7": 7,"chr8": 8, "chr9": 9, "chr10": 10, "chr11": 11, "chr12": 12, "chr13": 13, "chr14":14,"chr15":15,"chr16": 16, "chr17": 17, "chr18": 18, "chr19": 19,"chr20":20,"chr21":21,"chr22":22,"chrX": 23,"chrY": 24,"chrM": 25}

BINDIR="../.."
APPSDIR=os.path.join(BINDIR + "/applications/")
SAMTOOLS=os.path.join(APPSDIR + "/samtools/samtools")
BWA=os.path.join(APPSDIR + "/bwa-mem2/bwa-mem2")
MINIMAP2=os.path.join(APPSDIR + "/mm2-fast/minimap2")

def populate_yaml(args):
    with open(args['config'], 'r') as f:
        data = yaml.safe_load(f)
    #print(data)
    for keys in ['dataset', 'fqprocess', 'bwa', 'mm2']:
        arg = data[keys]  #['arguments']
        for k,v in arg.items():
            #print('k: ', k, args[k], ' keys: ', keys)
            if v != "None" and v != "":
                #print(k,v)
                if k == "params":
                    if keys == 'bwa' and args['read_type'] == 'short':
                        args[k] = v
                    if keys == 'mm2' and args['read_type'] == 'long':
                        args[k] = v
                    #else:
                    #    assert True, 'Error: Incorrect parameter settings in config: read_type/params'
                else:
                    args[k] = v

    return args

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
    #print('k: ', k, " d: ", d)
        
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
    #print('unlink the pipe...')

    
def fastq_reader_real( comm, cpus, fname, outpipe, outdir, last=False ):
    ## TBW
    pass
    f = open(fname, "r")
    f.seek(0, 2) # go to end
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
        print("[Info] Time for pragzip index: "+str(t1-t0))
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
chromo_dict={}
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
    sz = 0
    while ndone < nranks:
        bf_sz = 100000000
        #bf_sz2 = 100 #100000
        #req1 = comm.irecv(bf_sz, tag=0)
        #stat1 = MPI.Status()
        #vl = req1.wait(stat1)
        #print(f'[{rank}] recv size vl: ', stat1.Get_count(MPI.CHAR))
        ##print(f'[{rank}] Error vl: ', MPI.Get_error_string(stat1.Get_error()))
        #req = comm.irecv(bf_sz2, tag=1)
        #stat = MPI.Status()
        #print(f'[{rank}] status: ', stat)
        #print(f'[{rank}] req', req)
        #b = req.wait(stat)
        #print(f'[{rank}] #size: ', stat.Get_count(MPI.CHAR))
        #req = comm.irecv()
        #b, vl = req.wait()
        req = comm.irecv(bf_sz)
        b, vl = req.wait()
        #recv_rank = b
        #print(f" [{rank}] after wait ", ndone, "/", nranks, ' recvd from: ', b)
        b = b//nranks
        if len(vl)>0 and  vl[0]=="done":
            ndone+=1
            #print(f"[{rank}] recvd done from rank ", recv_rank, sz)
            #print ("sort_thr "+str(rank)+", got done:"+str(ndone))
        else:
            fl[b].writelines(vl)
            #pass
    [f.close() for f in fl]
    t1 = time.time()
    #print(f"[{rank}] Done sort_thr "+str(rank)+" time: "+str(t1-t0))


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
            #assert len(send_bufs[b]) < 10000000
            req = comm.isend( (b, send_bufs[b]),r)
            req.wait()
            #req.test()
            #req1 = comm.isend(send_bufs[b],r, tag=0)
            ##req1 = comm.isend(["wasim"],r)
            #stat1 = MPI.Status()
            #req1.wait(stat1)
            #print(f"[{rank}] wait", stat1.Get_count(MPI.CHAR), len(send_bufs[b]))
            #assert stat1.Get_count(MPI.CHAR) < 100000000
            #
            #req = comm.isend(b,r, tag=1)
            #stat = MPI.Status()
            #req.wait(stat)
            #print(f"[{rank}] free", stat.Get_count(MPI.CHAR))

            #if not req.Test():
            #    sent_msgs.append( (req, send_bufs[b]) )  # if send is not yet complete, save this request and data
            #    print("saved request",rank,r)
            #req.wait()
            send_bufs[b].clear()
            #send_bufs[b] = []
            #if force: print ("send forced on rank "+str(rank)+" to "+str(r))
            return req
        #sent_msgs = [ rqt for rqt in sent_msgs if not rqt[0].Test() ]
    def read_wrap(f):
        return f.readline()

    t0 = time.time()
    #seq_start['*'] = 0    # put all unmatched reads at beginning
    cumlen = 0
    global keep
    global chromo_dict
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
                    #if sn in mydict.keys():
                    if sn in chromo_dict.keys():
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

        #print(f"[{rank}] done part 1 ")            
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
        if r!=rank:
            send_wrap("done", r, force=True)
            #print(f"[{rank}] sending done to rank: ", r)
    send_wrap("done", rank, force=True)  # send to self last to avoid race condition with allreduce
    t1 = time.time()
    #cnt=0
    #while cnt<4:
    #    for rq in rt:
    #        flag, status = rq.test()
    #        if flag:
    #            cnt+=1
    #            print(status, cnt, flush=True)
    #print(f"{rank} check", flush=True)
    #total = comm2.allreduce(i)  ## put it after confirming that all the send and recv done
    #total = comm.allreduce(i)
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
    parser.add_argument('--tempdir',default="",help="Intermediate data directory")
    parser.add_argument('--refdir',default="",help="Reference genome directory")
    parser.add_argument('--read1',default="", nargs='+',help="name of r1 files")
    parser.add_argument('--read2',default="", nargs='+',help="name of r2 files")
    parser.add_argument('--read3',default="", nargs='+',help="name of r3 files (for fqprocess) seperated by spaces")
    parser.add_argument('--readi1',default="", nargs='+',help="name of i1 files (for fqprocess) seperated by spaces")
    
    parser.add_argument('--prefix',default="", help="prefix for processed R1 and R3 files for bwa-mem2")
    parser.add_argument('--suffix',default="", help="suffix for processed R1 and R3 files for bwa-mem2")
   
    parser.add_argument('--whitelist',default="whitelist.txt",help="10x whitelist file")
    parser.add_argument('--read_structure',default="16C",help="read structure")
    parser.add_argument('--barcode_orientation',default="FIRST_BP_RC",help="barcode orientation")
    parser.add_argument('--sample_id',default="",help="sample id")
    parser.add_argument('--output_format',default="",help="output_format")
    parser.add_argument('--output',help="Output data directory")
    parser.add_argument('--mode', default='pragzip', help="flatmode/fqprocessonly/multifq/pragzip. flatmode is just bwa w/o sort.")
    parser.add_argument('--params', default='', help="parameter string to bwa-mem2 barring threads paramter")
    parser.add_argument("-i", "--index", help="name of index file")
    parser.add_argument("-p", "--outfile", help="prefix for read files")
    #parser.add_argument("-r", "--preads", nargs='+',help="name of reads file seperated by space")
    parser.add_argument("-c", "--cpus",default=1,help="Number of cpus. default=1")
    parser.add_argument("-t", "--threads",default=1,help="Number of threads used in samtool operations. default=1")
    parser.add_argument("-b", "--bam_size",default=9, help="bam file size in GB")
    parser.add_argument('-in', '--rindex',action='store_true',help="It will index reference genome for bwa-mem2. If it is already done offline then don't use this flag.")
    parser.add_argument('-dindex',action='store_true',help="It will create .fai index. If it is done offline then disable this.")
    #parser.add_argument('--container_tool',default="docker",help="Container tool used in pipeline : Docker/Podman")
    #parser.add_argument('--shards',default=1,help="Number of shards for deepvariant")
    parser.add_argument('-pr', '--profile',action='store_true',help="Use profiling")
    parser.add_argument('--keep_unmapped',action='store_true',help="Keep Unmapped entries at the end of sam file.")
    #parser.add_argument('--se_mode',action='store_true',help="Single End (SE) reads only. Paired End (PE) mode is default")
    parser.add_argument('--read_type',default="short", help="(short/long): bwa-mem2 with short reads, mm2-fast with long reads")
    parser.add_argument('-y', "--config", default="", help="config yaml file for args")
    ## Arg values in the provided yaml file override any command-line or default args values
    args = vars(parser.parse_args())

    #print('args: ', args)
    if args["config"] != "":
        args = populate_yaml(args)

    #print('args: ', args)
    global chromo_dict
    #ncpus = int(cpus)
    #start0 = time.time()
    
    #se_mode=args["se_mode"]
    read_type=args["read_type"]
        
    ifile=args["index"]
    params=args["params"]

    params=params.replace("+","-")
    read1 = rfile1 = args["read1"]
    read2 = rfile2 = args["read2"]
    read3 = args["read3"]
    readi1 = args["readi1"]
    se_mode = False
    if read2 == "" or read2 == "None":
        se_mode = True
        
    if read_type == "long":
        #assert se_mode == True, "for long reads se_mode should be enabled by the code"
        se_mode = True   ## for long reads, se_mode=True always
    
    whitelist=args["whitelist"]
    read_structure=args["read_structure"]
    barcode_orientation=args["barcode_orientation"]
    bam_size=args["bam_size"]
    outfile=args["outfile"]
    mode=args["mode"]
    cpus=args["cpus"]
    threads=args["threads"]    ## read prefix for R1, I1, R2 files
    rindex=args["rindex"]
    folder=args["input"] + "/"
    output=args["output"] + "/"
    tempdir=args["tempdir"]
    if tempdir=="": tempdir=output
    refdir=args["refdir"]
    #if refdir=="": refdir=folder
    prof=args["profile"]
    global keep
    keep=args["keep_unmapped"]
    dindex=args["dindex"]
    
    sample_id=args['sample_id']
    if sample_id == "": sample_id = output
    output_format = args["output_format"]
    
    prefix=args["prefix"]  ## prefix for mutlifq2sortedbam mode reading 'fqprocess' processed fastq files
    suffix=args["suffix"]  ## suffix for mutlifq2sortedbam mode reading 'fqprocess' processed fastq files
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()

    #print(refdir + ifile)
    #global bins_per_rank
    #bins_per_rank = max(4,bins_per_rank//nranks)
    global ncpus
    ncpus = int(cpus)
    
    if rank == 0:
        #print("\n")
        if se_mode:
            print("[Info] Running in SE reads only mode!!!")
        else: print("[Info] Running in PE reads mode!!!")
        if read_type == "short":
            print("[Info] BWA-MEM2 parameters: ", params)
        if read_type == "long":
            print("[Info] mm2-fast parameters: ", params)
    

    def faidx(refdir, ifile):
        tic=time.time()
        #print(f'{SAMTOOLS} faidx '+os.path.join(refdir,ifile) + ' > samlog.txt')
        a=run(f'{SAMTOOLS} faidx '+os.path.join(refdir,ifile) + ' > ' + output + 'logs/samlog.txt', capture_output=True, shell=True)
        assert a.returncode==0,"samtools reference index creation failed"
        end=time.time()
        print("\n[Info] Reference to .fai index creation time",end-tic)
        
    #start0 = time.time()
    # Generate data (reads) index file(s)

    flg=0
    if rank == 0:
        if dindex == 'True':
            #print("############# dindex: ", dindex)
            g = refdir + ifile
            print(g.split(".")[-1])
            if g.split(".")[-1] != "gz":
                faidx(refdir, ifile)
            else:
                print("Error: genome file should not be gzip for indexing.\nExiting...")
                flg=1

    comm.barrier()
    flg = comm.bcast(flg, root=0)
    if flg:
        os.sys.exit(1)
        #comm.MPI_Finalize()
        

    if rank==0:
        #yappi.set_clock_type("wall")
        #if prof: yappi.start()
        #file_size = os.path.getsize(folder+rfile1)
        #print("\nSize of FASTQ file:",file_size)
        if rindex == 'True' :
            print("[Info] Indexing Starts", flush=True)
            begin = time.time()
            a=run(f'{BWA} index '+ refdir + ifile + ' > ' + output + '/logs/bwaindexlog.txt', capture_output=True,shell=True)
            end=time.time()
            #file_size = os.path.getsize(folder+rfile1)
            print("\n[Info] Bwa Index creation time:",end-begin)
            #print("\n[Info] Size of FASTQ file:",file_size)

    comm.barrier()             
    #################################################################
    if mode in ['flatmode', 'pragzip', 'multifq']:
        ## bwa-mem2 index path check, every ranks checks the access.
        #g = refdir + ifile
        #print(g.split(".")[-1])
        #if g.split(".")[-1] != "gz":
        #  if rank == 0:
        #    faidx(refdir, ifile)
        #else:
        #    print("Error: genome file should not be gzip. Exiting...")
        #    os.sys.exit(1)
        #comm.barrier()
        
        #print("Genome: ", refdir + ifile)
        assert os.path.exists(refdir + ifile) == True, "missing bwa-mem2 genome"
        if not read_type == "long":
            flg = 0
            #assert os.path.exists(refdir + ifile + ".bwt.2bit.64") == True, "missing bwa-mem2 index files"
            #assert os.path.exists(refdir + ifile + ".0123") == True, "missing bwa-mem2 index files"
            #assert os.path.exists(refdir + ifile + ".amb") == True, "missing bwa-mem2 index files"
            #assert os.path.exists(refdir + ifile + ".ann") == True, "missing bwa-mem2 index files"
            #assert os.path.exists(refdir + ifile + ".pac") == True, "missing bwa-mem2 index files"
            if os.path.exists(refdir + ifile + ".bwt.2bit.64") == False: flg=1
            if os.path.exists(refdir + ifile + ".0123") == False: flg = 1
            if os.path.exists(refdir + ifile + ".amb") == False: flg=1
            if os.path.exists(refdir + ifile + ".ann") == False: flg=1
            if os.path.exists(refdir + ifile + ".pac") == False: flg=1
            if flg and rank == 0:
                print("Error: Some BWA reference index files missing...Exiting")
                print("Error: You can enable rindex option in the config file to create index on-the-fly.")
                print("Error: Or, you can manually build and place the index in ", refdir)
            if flg:
                os.sys.exit(1)
            

    if mode in ['flatmode', 'pragzip']:
        #print("Read1: ", folder + rfile1)
        assert os.path.exists(folder+rfile1) == True, "missing input read1 files"
        if not se_mode:
            #print("Read2: ", folder + rfile2)
            assert os.path.exists(folder+rfile2) == True, "missing input read2 files"
        # output paths
        #print("Output dir: ", folder)
        assert os.path.exists(output) == True, "output path does not exist"

    if mode in ['pragzip', 'flatmode']:
        ## bwa-mem2 index path check, every ranks checks the access.
        g = rfile1
        if g.split(".")[-1] != "gz":
            print("Error: reads file1 not in gzip format. Exiting...")
            os.sys.exit(1)
            
        if not se_mode:
            g = rfile2
            if g.split(".")[-1] != "gz":
                print("Error: reads file2 not in gzip format. Exiting...")
                os.sys.exit(1)
            
        comm.barrier()

    # chromo_dict information added
    #chromo_file=os.path.join(refdir,ifile)+".fai"
    #if os.path.isfile(chromo_file):
    #    f = open(chromo_file, "r")
    #else :
    #    print(f'Error: File "{chromo_file}" not found')
    #    exit()
    #    
    #lines=f.readlines()
    #f.close()
    #for line in lines[0:25]:
    #    chromo=line.split("\t")[0]
    #    chromo_dict[chromo]= 1
        
    if rank == 0:
        fo = os.path.join(output, "logs")
        if not os.path.exists(fo):
            try:
                os.mkdir(fo)
            except:
                print("Error: Unable to create logs folder inside the ouptut folder")
                os.sys.exit(1)
        else:
            print('[Info] output logs folder exits, will override the log files')
            
    if mode == 'flatmode':
        '''
        Goal: Takes the reads (gzip) SE file or PE files as input and 
              produces unsorted SAM files as output
        Input: 
        - bwa-mem2 index files
        - SE reads file or PE reads files

        Output:
        - Outputs multiple unsroted sam files prefixed with "aln" in the output folder

        Attributes:
        - It produces mutlitple, equal to the number of ranks, SAM files as the ouput
        - It does not sort or merge the ouptut SAM files
        - It supports SE reads (1 read files) or PE reads (2 read files) as input
        - It needs reads files in gzip format
        '''
        
        #if rank == 0:
        #    print("dist bwa-mem2 starts in flatmode")
        if rank == 0 and read_type == "short":
            print("[Info] bwa-mem2 starts in flatmode")
        if rank == 0 and read_type == "long":
            print("[Info] mm2-fast starts in flatmode")
            
        # Execute bwamem2 -- may include sort, merge depending on mode
        begin0 = time.time()
        fn1 = pragzip_reader( comm, int(cpus), folder+rfile1, output )
        if not se_mode:
            fn2 = pragzip_reader( comm, int(cpus), folder+rfile2, output, last=True )
        fn3 = os.path.join(output, 'aln' + str("%05d"%rank) + ".sam")
        begin1 = time.time()
        #a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3,capture_output=True, shell=True)
        #bwastr = '{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'/bwalog' + str(rank) + '.txt'
        #print(bwastr)
        if se_mode:
            if read_type == "long":
                ## s='minimap2 ' + params  + ' -t '+cpus+' '+refdir+ifile+' '+fn1+' '+' >     '+fn3 + '  2> ' + output +'logs/mm2log' + str(rank) + '.txt'
                #print('mm2 cmd: ', s, flush=True)
                a=run(f'{MINIMAP2} ' + params + ' -t '+cpus+' '+refdir+ifile+' '+fn1+' '+' > '+fn3 + '  2> ' + output +'logs/mm2log' + str(rank) + '.txt',capture_output=True, shell=True)
            else:
                a=run(f'{BWA} mem ' + params + ' -t '+cpus+' '+refdir+ifile+' '+fn1+' '+' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)

        else:                
            a=run(f'{BWA} mem ' + params + ' -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)
            
        assert a.returncode == 0
        
        end1b=time.time()
        #thr.join()
        comm.barrier()
        end1=time.time()

        if rank==0:
            print("\n[Info] FASTQ to SAM time (flatmode):",end1-begin1)
            print("   (includes wait time:",end1-end1b,")")
            print("[Info] Output SAM files are present in the ouptut folder")
            #print("\nsam to sort-bam starts")
            #begin2=time.time()

        return

    elif mode == 'fqprocessonly':
        ## fastq preprocess and split
        if rank == 0:
            print("#############################################")
            print("[Info] Whitelist: ", whitelist)
            print("[Info] read_structure: ", read_structure)
            print("[Info] barcode_orientation: ", barcode_orientation)
            print("[Info] bam_size: ", bam_size)
            print("[Info] sample_id: ", sample_id)
            print("[Info] output_format: ", output_format)
            print("#############################################")
            print("[Info] Starting fqprocess....", flush=True)

            tic = time.time()
            #rl1, rl2, rl3 = read1.split(), read2.split(), read3.split()
            rl1, rl2, rl3, ri1 = read1, read2, read3, readi1
            barcode_index = rl2[0]
            #cmd='zcat ' + folder/barcode_index '| sed -n "2~4p" | shuf -n 1000 > ' + output + '/downsample.fq'
            bash_command = '''zcat '''+  folder + "/" + barcode_index + \
            ''' | sed -n '2~4p' | shuf -n 1000 > downsample.fq'''
            a = run(bash_command,shell=True, check=True, executable='/bin/bash')
            assert a.returncode == 0

            # dynamic-barcode-orientation downsample.fq  whitelist best_match.txt
            try:
                cmd='python3 warp-tools/tools/scripts/dynamic-barcode-orientation.py downsample.fq ' \
                    + whitelist + ' best_match.txt'
                a = run(cmd, shell=True, check=True, executable='/bin/bash')
                assert a.returncode == 0
            except:
                print("[Info] Exception in dynamic-barcode-orientation!!")

            # Read the contents of the file into a variable
            with open('best_match.txt', 'r') as file:
                barcode_choice = file.read().strip()

            #r1, r2, r3="--R1 ", "--R2 ", "--R3 "
            r1, r2, r3, i1 ="", "", "", ""
            for r in range(len(rl1)):
                r1+="--R1 " + folder + "/" + rl1[r] + ' '
                r2+="--R2 " + folder + "/" + rl2[r] + ' '
                r3+="--R3 " + folder + "/" + rl3[r] + ' '

            for r in range(len(ri1)):
                i1+="--I1 " + folder + "/" + ri1[r] + ' '

            ## print(r1)
            ## print(r2)
            ## print(r3)
            cmd=f'warp-tools/tools/fastqpreprocessing/bin/fastqprocess --verbose --bam-size ' \
                +  str(bam_size) + \
                ' --sample-id ' + sample_id + ' ' + r1 + ' ' + r2 + ' ' + r3 + \
                ' --output-format ' + output_format +  '  --barcode-orientation ' + barcode_choice + \
                ' --read-structure ' + read_structure + '  --white-list ' + whitelist

            print('[Info] fastqprocess cmd: ', cmd)
            a=run(cmd, shell=True, check=True, executable='/bin/bash')
            assert a.returncode == 0
            ## Create #output files equal to the number of ranks
            ## create output file names as rprefix + str(rank) + '_1/2.fastq.gz'
            toc = time.time()
            print("[Info] Time for fq processing: ", toc - tic)

        # Preindex reference genome if requested
        comm.barrier()
        return

    elif mode == "multifq":
        if rank==0:
            print("\n[Info] bwa-mem2 starts..")

        # Execute bwamem2 -- may include sort, merge depending on mode
        begin0 = time.time()
        if rank == 0:
            for r in range(nranks):
                #fn1 = "fastq_R1_" + str(r) + ".fastq.gz"
                #fn2 = "fastq_R3_" + str(r) + ".fastq.gz"
                #multiome-practice-may15_arcgtf_0.R1.trimmed_adapters.fastq.gz
                fn1 = os.path.join(folder, prefix + "_" + str(r) + ".R1." + suffix)
                fn2 =  os.path.join(folder, prefix + "_" + str(r) + ".R3." + suffix)

                if os.path.isfile(fn1) == False or os.path.isfile(fn2) == False:
                    print(f"[Info] Error: Number of files fastq files ({r}) < number of ranks ({nranks})")
                    print(f"!!! Fastq file(s) are not available for processing by rank {r} aborting..\n\n")
                    #sys.exit(1)
                    error_code = 1
                    comm.Abort(error_code)

        comm.barrier()
        #fn1 = "fastq_R1_" + str(rank) + ".fastq.gz"
        #fn2 = "fastq_R3_" + str(rank) + ".fastq.gz"
        fn1 = os.path.join(folder, prefix + "_" + str(rank) + ".R1." + suffix)
        fn2 = os.path.join(folder, prefix + "_" + str(rank) + ".R3." + suffix)

        #fn1 = folder + "/fastq_R1_" + str(rank) + ".fastq.gz"
        #fn2 = folder + "/fastq_R3_" + str(rank) + ".fastq.gz"
        if rank == 0:
            print("[Info] Input files: ")
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
            a=run(f'{BWA} mem ' + params + ' -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)
            assert a.returncode == 0
        else:
            print(f"[Info] {rank} No input file for me")
        end1b=time.time()

        thr.join()
        comm.barrier()
        end1=time.time()

        if rank==0:
            print("[Info] \nFASTQ to SAM time (fqprocess):",end1-begin1)
            print("   (includes wait time:",end1 - end1b,")")

            print("\n[Info] sam to sort-bam starts")
            gin2=time.time()

            
    elif mode == 'pragzip':
        # Preindex refernce genome if requested
        '''
        Goal: Takes the reads (gzip) SE file or PE files as input and 
              produces sorted BAM files as output
        Input: 
        - bwa-mem2 index files
        - SE reads file or PE reads files

        Output:
        - Outputs BAM file of name

        Attributes:
        - It produces one sorted bam file
        - It supports SE reads (1 read files) or PE reads (2 read files) as input
        - It needs reads files in gzip format, unzip files support WIP. Let me know if you need it.
        '''

        if rank == 0 and read_type == "short":
            print("[Info] bwa-mem2 starts")
        if rank == 0 and read_type == "long":
            print("[Info] mm2-fast starts")

        # Execute bwamem2 -- may include sort, merge depending on mode
        begin0 = time.time()
        #try:
        fn1 = pragzip_reader( comm, int(cpus), folder+rfile1, output )
        if not se_mode:
            fn2 = pragzip_reader( comm, int(cpus), folder+rfile2, output, last=True )
        fn3, thr = sam_writer( comm, output+'/aln' )
        #fn3 = output + "/aln" + str(rank) + ".sam"
        begin1 = time.time()
        if se_mode:
            if read_type == "long":
                a=run(f'{MINIMAP2} ' + params + ' -t '+cpus+' '+refdir+ifile+' '+fn1+' '+' > '+fn3 + '  2> ' + output +'logs/mm2log' + str(rank) + '.txt',capture_output=True, shell=True)
            else:                
                a=run(f'{BWA} mem ' + params +' -t '+cpus+' '+refdir+ifile+' '+fn1+' '+' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)
        else:
            a=run(f'{BWA} mem ' + params +' -t '+cpus+' '+refdir+ifile+' '+fn1+' '+fn2+' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)
            
        assert a.returncode == 0
        end1b=time.time()
        thr.join()
        comm.barrier()
        end1=time.time()
            #print('end: ', end, flush=True)
        #except Exception as e:
        #    print(f"An error occurred: {e}")
        #    print("Exiting...")
        #    print("Please check bwa logs in the output folder for more details.")
        #    os.sys.exit(1)

        if rank==0:
            print("\n[Info] FASTQ to SAM time:",end1-begin1)
            print("   (includes wait time:",end1-end1b,")")

            print("\n[Info] sam to sort-bam starts")
            begin2 = time.time()

    else:
        print("Error: Incorrect mode specified!\n   please specify: pragzip/flatmode/multifq/fqprocessonly")
        os.sys.exit(1)

    # Finish sort, merge, convert to bam depending on mode
    cmd=""
    for i in range(bins_per_rank):
        binstr = '%05d'%(nranks*i+rank)
        cmd+=f'{SAMTOOLS} sort --threads '+threads+' -T '+tempdir+'/aln'+binstr+'.sorted -o '+ output +'/aln'+binstr+'.bam '+ output+'/aln'+binstr+'.sam;'
        if i%20==0:
            a=run(cmd,capture_output=True,shell=True)
            cmd=""
            
    if not cmd=="": a=run(cmd,capture_output=True,shell=True)
    comm.barrier()

    if rank==0:
        end2=time.time()
        print("[Info] SAM to sort-BAM time:",end2-begin2)

    ## concat bams
    if rank == 0:
        tic = time.time()
        bf = []
        print('[Info] Concating the bam files...')
        for b in range(bins_per_rank):
            for r in range(nranks):
                binstr = '%05d'%(nranks*b + r)
                bf.append(output+'/aln'+binstr+'.bam')

        infstr = bf[0]
        for i in range(1, len(bf)):
            infstr = infstr + " " + bf[i]
        if outfile == None:
            outfile = "final"

        cmd+=f'{SAMTOOLS} cat -o ' + os.path.join(output, outfile) + '.sorted.bam ' + infstr
        #print("merge cmd: ", cmd, flush=True)
        a=run(cmd,capture_output=True,shell=True)
        assert a.returncode == 0
        print("[Info] Concat done, time taken: ", time.time() - tic)
        
    if rank == 0:
        print("[Info] Cleaning up...")
        tic = time.time()
        os.system('rm ' + output + "*.sam")
        os.system('rm ' + output + "aln*.bam")
        os.system('rm ' + output + "*.idx")
        toc = time.time()        
        print("[Info] Clean up done, time taken: ", toc -tic)

        
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
