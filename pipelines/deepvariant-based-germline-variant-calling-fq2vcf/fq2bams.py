#*************************************************************************************
#                           The MIT License
#
#   Intel OpenOmics - fq2vcf pipeline
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
#Authors: Manasi Tiwari <manasi.tiwari@intel.com>;  Vasimuddin Md <vasimuddin.md@intel.com>; Babu Pillai <padmanabhan.s.pillai@intel.com>; Sanchit Misra <sanchit.misra@intel.com>
#*****************************************************************************************/

from subprocess import Popen, PIPE, run
import subprocess
import json, time, os, sys
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
import pickle

BINDIR="../.."
APPSDIR=os.path.join(BINDIR + "/applications/")
SAMTOOLS=os.path.join(APPSDIR + "/samtools/samtools")
BWA=os.path.join(APPSDIR + "/bwa-mem2/bwa-mem2")
BWASSE=os.path.join(APPSDIR + "/bwa-mem2/bwa-mem2.sse42")

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
    idxname = os.path.join(outdir,fname.split('/')[-1]+".idx")
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
chromo_dict={}
# NOTE:  Code relies on insert order of dictionary keys -- so needs python>=3.7

# Calculate bins
def calculate_bins(nranks):
    global bins
    global bin_region
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
    #a = run('ls -las fname+("%05d.sam"%(0*nranks+rank))', capture_output=False, shell=True)
    #s =  fname+("%05d.sam"%(0*nranks+rank))
    #print(s)

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
                    #print(sn)
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
    #total = comm2.allreduce(i)
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


def allexit(comm, flg):
    comm.barrier()        
    flg = comm.bcast(flg, root=0)
    if flg: os.sys.exit(1)


def faidx(refdir, ifile, output):
    g = refdir + ifile
    #print(g.split(".")[-1])
    if g.split(".")[-1] == "gz":
        print("Error: genome file should not be gzip for indexing.\nExiting...")
        return 1
    
    tic=time.time()
    a=run(f'{SAMTOOLS} faidx '+os.path.join(refdir,ifile) + ' > ' + \
          output + 'logs/samfaidxlog.txt', capture_output=True, shell=True)
    assert a.returncode==0,"samtools reference index creation failed"
    end=time.time()
    print("\n[Info] Reference to .fai index creation time",end-tic)
    return 0


def input_check(rank, comm, refdir, inputdir, output, se_mode, ifile, rfile1, rfile2):
    assert os.path.exists(refdir + ifile) == True, "missing bwa-mem2 genome"
    flg = 0
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

    assert os.path.exists(inputdir+rfile1) == True, "missing input read1 files"
    if not se_mode:
        assert os.path.exists(inputdir+rfile2) == True, "missing input read2 files"
    # output paths
    assert os.path.exists(output) == True, "output path does not exist"

    flg = 0
    if rank == 0:
        g = rfile1
        if g.split(".")[-1] != "gz":
            print("Error: reads file1 not in gzip format. Exiting...")
            #os.sys.exit(1)
            flg = 1
            
        if not se_mode:
            g = rfile2
            if g.split(".")[-1] != "gz":
                print("Error: reads file2 not in gzip format. Exiting...")
                #os.sys.exit(1)
                flg = 1
                    
    allexit(comm, flg)  ## all ranks exit if failure in rank 0 above       

def create_folder(fo):
    flg = 0
    if not os.path.exists(fo):
        try:
            os.mkdir(fo)
        except:
            print("Error: Unable to create logs folder inside the ouptut folder")
            flg = 1
    else: print('[Info] output logs folder exits, will override the log files')
    return flg


def rundown(args):
    ifile = args["refindex"]
    rfile1 = args["read1"]
    se_mode = True
    rfile2 = ""
    if args["read2"] != "" and args["read2"] != "None":
        rfile2 = args["read2"]
        se_mode = False

    params=args["params"]
    #params = ""
    #if args["params"] != "" and args["params"] != "None":
    #    params = " -R " +  "\"" + args["params"] + "\""
    
    cpus = args["cpus"]
    threads = args["threads"]
    rindex = args["rindex"]
    dindex = args["dindex"]
    #nproc = args["shards"]
    inputdir = args["input"] + "/"
    output = args["output"] + "/"
    tempdir = args["tempdir"]
    if tempdir == "" or tempdir == "None": tempdir = output
    else: tempdir = tempdir + "/"
    refdir = args["refdir"] + "/"
    #if tempdir=="": tempdir=output
    #refdir=args["refdir"]
    #if refdir=="": refdir=inputdir
    prof=args["profile"]
    global keep
    #keep=args["keep_unmapped"]
    if not args["not_keep_unmapped"]: keep=True
    keep_sam=args["keep_intermediate_sam"]
    container_tool=args["container_tool"]
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()

    #global bins_per_rank
    #bins_per_rank = max(4,bins_per_rank//nranks)
    global ncpus
    global chromo_dict
    ncpus = int(cpus)

    #print('tempdir: ', tempdir)
    if rank == 0:
        if se_mode:
            print("[Info] Running in SE reads only mode!!!")
        else: print("[Info] Running in PE reads mode!!!")

    ########################################################################            
    #start0 = time.time()
    flg = 0
    if rank == 0:
        fo = os.path.join(output, "logs")
        flg = create_folder(fo)
    allexit(comm, flg)        
    ########################################################################    
    
    start0 = time.time()
    # chromo_dict information added
    ## chromo_file=os.path.join(refdir,ifile)+".fai"
    ## if os.path.isfile(chromo_file):
    ##     f = open(chromo_file, "r")
    ## else :
    ##     print(f'File "{chromo_file}" not found')
    ##     exit()
    ## lines=f.readlines()
    ## f.close()
    ## for line in lines[0:25]:
    ##     chromo=line.split("\t")[0]
    ##     chromo_dict[chromo]= 1
    
    # Generate data (reads) index file(s)
    flg=0
    if rank == 0 and dindex == 'True':
        flg = faidx(refdir, ifile)
        if flg:
            print("[Info] Exiting due to a failure in creating .fai file.")            
    allexit(comm, flg)
        
    # Preindex refernce genome if requested
    flg = 0
    if rank==0 and rindex == 'True':
        print("[Info] Indexing Starts", flush=True)
        begin = time.time()
        flg = run(f'{BWA} index '+ refdir + ifile + ' > ' + output + \
                  '/logs/bwaindexlog.txt', capture_output=True,shell=True)
        flg = flg.returncode
        end=time.time()
        print("\n[Info] Bwa Index creation time:", end - begin)
        if flg:
            print("[Info] BWA-MEM2 index creation failed. Check logs.")
    allexit(comm, flg)                    

    ##############################################################################
    input_check(rank, comm, refdir, inputdir, output, se_mode, ifile, rfile1, rfile2)
    ##############################################################################
    
    # chromo_dict information added
    if keep == False: ## and rank == 0:
        chromo_file = os.path.join(refdir,ifile)+".fai"
        if os.path.isfile(chromo_file):
            f = open(chromo_file, "r")
            lines=f.readlines()
            f.close()
            for line in lines[0:25]:
                chromo=line.split("\t")[0]
                chromo_dict[chromo]= 1
            #print(f'chromo_dict: {chromo_dict}')
        else:
            if rank == 0:
                print(f'[Info] File "{chromo_file}" not found, resetting keep=True')
            #exit()
            keep = True                


    comm.barrier()                
    #keep = comm.bcast(keep, root=0)            
    #############################################################################
    if rank == 0:
        print("[Info] bwa-mem2 starts")
    
    begin0 = time.time()
    fn1 = pragzip_reader( comm, int(cpus), os.path.join(inputdir,rfile1), output )
    if not se_mode:    
        fn2 = pragzip_reader( comm, int(cpus), os.path.join(inputdir,rfile2), output, last=True )
    fn3, thr = sam_writer( comm, os.path.join(tempdir,'aln') )
    #print("tempdir:" , tempdir, os.path.join(tempdir,'aln'))
    begin1 = time.time()
    #a=run(f'{BINDIR}/applications/bwa-mem2/bwa-mem2 mem -t '+cpus+' '+os.path.join(refdir,ifile)+' '+fn1+' '+fn2+' > '+fn3,capture_output=True, shell=True)
    if se_mode:
        if args['simd'] == 'sse':
            cmd = f'{BWASSE} mem ' + params + ' -t '+cpus+' '+ refdir+ifile+' '+fn1+' '+' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) + '.txt'
        else:
            cmd = f'{BWA} mem ' + params + ' -t '+cpus+' '+ refdir+ifile+' '+fn1+' '+' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) + '.txt'
            
        #print(cmd)
        a=run(cmd, capture_output=True, shell=True)
    else:
        if args['simd'] == 'sse':
            a=run(f'{BWASSE} mem ' + params + ' -t '+cpus+' '+refdir+ifile+' '+ \
              fn1+' '+fn2+' > '+fn3 + '  2> ' + output + \
                  'logs/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)
            
        else:                
            a=run(f'{BWA} mem ' + params + ' -t '+cpus+' '+refdir+ifile+' '+ \
                  fn1+' '+fn2+' > '+fn3 + '  2> ' + output + \
                  'logs/bwalog' + str(rank) + '.txt',capture_output=True, shell=True)
    
    assert a.returncode==0,"bwa-mem2 Run Failed"
    end1b=time.time()
    thr.join()
    comm.barrier()
    end1=time.time()

    if rank==0:
        print("\n[Info] FASTQ to SAM time:",end1-begin1)
        print("   (includes wait time:",end1-end1b,")")        
        print("\n[Info] SAM to sort-BAM starts")
        begin2 = time.time()

    ##############################################################################
    # Finish sort, merge, convert to bam depending on mode
    cmd=""
    for i in range(bins_per_rank):
        binstr = '%05d'%(nranks*i+rank)
        #cmd+=f'{BINDIR}/applications/samtools/samtools sort --threads '+threads+' -T '+os.path.join(tempdir,'aln'+binstr+'.sorted') + ' -o '+ os.path.join(tempdir,'aln'+binstr+'.bam')+' '+ os.path.join(tempdir,'aln'+binstr+'.sam')
        cmd+=f'{SAMTOOLS} sort --threads '+threads+' -T '+tempdir+ '/aln'+binstr+ \
            '.sorted -o '+ tempdir +'/aln'+binstr+'.bam '+ tempdir+'/aln'+binstr+'.sam' + " > " + output + "/logs/samsortlog" + binstr + ".txt;"        
        #print('samsort: ', cmd)
        if i%20==0:
            a = run(cmd, capture_output = True, shell = True)
            assert a.returncode==0,"samtools sort Failed"
            cmd=""
            
    if not cmd=="":
        a=run(cmd,capture_output=True,shell=True)
        assert a.returncode==0,"samtools sort Failed"
    comm.barrier()

    
    if rank==0:
        end2=time.time()
        #print("SAM to sort-BAM time:",end2-begin2)
        print("[Info] SAM to sort-BAM time:",end2-begin2)        
        cmd=""

    ##############################################################################
    if not keep_sam:
        for i in range(bins_per_rank):
            binstr = '%05d'%(nranks*i+rank)
            cmd+=f'rm -rf '+os.path.join(tempdir,'aln'+binstr+'.sam')

        a = run(cmd, capture_output=True, shell=True)
        assert a.returncode==0,"[Info] Deletion of SAM file failed"

    ##############################################################################
    # Generate index file(s)
    begin3=time.time()
    cmd =""
    for i in range(bins_per_rank):
        fname = "aln%05d.bam"%(i*nranks+rank)
        #a=run(f'{BINDIR}/applications/samtools/samtools index -M -@ '+threads+' '+inputdir+fname,capture_output=True,shell=True)
        cmd += f'{SAMTOOLS} index -M -@ '+threads+' '+os.path.join(tempdir, fname)+';'
        if i%20 == 0:
            a = run(cmd,capture_output=True,shell=True)
            assert a.returncode==0, "[Info] samtools bam index creation failed"
            cmd=""
            
    if not cmd=="":
        a=run(cmd,capture_output=True,shell=True)
        assert a.returncode==0,"[Info] samtools bam index creation failed."
    comm.barrier()

    ##############################################################################
    if rank==0:
        end3=time.time()
        print("\n[Info] BAM to .bai index creation time",end3-begin3)
        #print("\n[Info] Starting Deepvariant execution...")
        begin5 = time.time()
        with open(os.path.join(tempdir,'bin_region.pkl'), 'wb') as f:
            global bin_region
            pickle.dump(bin_region, f)

    comm.barrier()
    if rank==0:
        end5=time.time()
        #print("[Info] Wrote intermediate files at ", tempdir)
        print("\n[Info] fq2bam runtime",end5-start0)
        #print("\n[Info] Starting Deepvariant execution...")
        idxname = os.path.join(output, rfile1 +".idx")
        #print("R1: ", idxname)
        os.system('rm ' + idxname)
        if not se_mode:
            idxname = os.path.join(output, rfile2 +".idx")
            #print("R2: ", idxname)
            os.system('rm ' + idxname)

    ##############################################################################
    return 0

#def main(argv):
def main(args):

    if args['buildindexonly']:
        ifile=args["refindex"]
        output = args["output"] + "/"
        refdir = args["refdir"] + "/"
        
        print("[Info] Indexing Starts", flush=True)
        flg = faidx(refdir, ifile, output)
        assert flg == 0, 'faidax failed.'

        begin = time.time()
        a=run(f'{BWA} index '+ os.path.join(refdir, ifile) + ' > ' + output + \
              '/logs/bwaindexlog.txt', capture_output=True, shell=True)
        
        end=time.time()
        print("\n[Info] Bwa Index creation time:",end-begin)
    else:
        rundown(args)
    
        
if __name__ == "__main__":
    #main(sys.argv[1:])
    args = json.loads(sys.argv[1])
    main(args)
