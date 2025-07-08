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
#!/usr/bin/env bash
import json, time, os, sys
from subprocess import Popen, PIPE, run
import subprocess
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
BWASSE=os.path.join(APPSDIR + "/bwa-mem2/bwa-mem2.sse42")
MINIMAP2=os.path.join(APPSDIR + "/mm2-fast/minimap2")
STAR=os.path.join(APPSDIR + "/STAR/source/STAR")

def populate_yaml(args):
    with open(args['config'], 'r') as f:
        data = yaml.safe_load(f)

    for keys in ['dataset', 'fqprocess', 'bwa', 'mm2']:
        arg = data[keys]  #['arguments']
        for k,v in arg.items():
            if v != "None" and v != "":
                if k == "params":
                    if keys == 'bwa' and args['read_type'] == 'short':
                        args[k] = v
                    if keys == 'mm2' and args['read_type'] == 'long':
                        args[k] = v
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

def fastq_reader(comm, cpus, fname, outdir, last=False):
    #outpipe = tempfile.mktemp()
    outpipe = tempfile.NamedTemporaryFile()
    temp=outpipe.name
    outpipe.close()
    #os.unlink(outpipe.name)
    outpipe=temp
    os.mkfifo(outpipe)
    comm2 = comm.Clone()
    threading.Thread(target=fastq_reader_real, args=(comm2, cpus, fname, outpipe, outdir, last)).start()
    return outpipe


def if_reader(comm, cpus, fname, outdir, last=False):
    print('fname: ', fname.split(".")[-1])
    if fname.split(".")[-1] == "gz":
        return pragzip_reader(comm, cpus, fname, outdir, last)
    else:
        return fastq_reader(comm, cpus, fname, outdir, last)   

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
        if args["read_type"] == "short":
            req = comm.irecv()
        else:
            bf_sz = 100000000
            req = comm.irecv(bf_sz)
        b, vl = req.wait()
        b = b//nranks
        if len(vl)>0 and  vl[0]=="done":
            ndone+=1
        else:
            fl[b].writelines(vl)

    [f.close() for f in fl]
    t1 = time.time()


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
            req = comm.isend( (b, send_bufs[b]),r)
            req.wait()
            send_bufs[b].clear()
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
                a = len(l[1].split(':')[0])
                sn = l[1][a+1:]
                ln = int(l[2].split(':')[1])
                if keep:
                    seq_start[sn] = cumlen
                    seq_len[sn] = ln
                    cumlen += ln
                    #print("In keep")
                else:
                    #if sn in mydict.keys():
                    if sn in chromo_dict.keys():
                        seq_start[sn] = cumlen
                        seq_len[sn] = ln
                        cumlen += ln
                    else :
                        seq_start[sn] = -1
            l = read_wrap(f1)
        # done reading headers
        if keep:
            seq_start['*'] = cumlen    # put all unmatched reads at end
        else:
            seq_start['*'] = -1

        calculate_bins(nranks)
        binstarts = [ b[0] for b in bins ]
        headers_done.release()   # allow output thread to begin

        # read remaining lines (aligned short reads), send to appropriate bin
        i = 0
        while l:
            x = l.split()
            seq = x[2]
            offset = int(x[3])
            if keep:
                key = seq_start[seq] + offset
                bn = bisect_wrap(binstarts, key) - 1
                _, r, b = bins[bn]
                send_wrap( l, r, b=bn )
            else :
                if seq_start[seq]!= -1:
                    key = seq_start[seq] + offset
                    temp2=time.time()
                    bn = bisect_wrap(binstarts, key) - 1
                    _, r, b = bins[bn]
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
    if flg:
        os.sys.exit(1)


def input_check(rank, comm, refdir, inputdir, output, se_mode, ifile, rfile1, rfile2, read_type):
    all_read_types = ['short', 'long', 'rna_short', 'meth']
    flg = 0
    if rank == 0:
        if all_read_types.count(read_type) != 1:
            print('[Error] Incorrect read_type: ', all_read_types)
            flg = 1
        #assert all_read_types.count(read_type) == 1, 'incorrect read_type provided.'
    allexit(comm, flg)        
    
    if read_type == 'short':
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


def alignment(p):
    outfile, read_type, params, refdir, ifile, fn1, fn2, fn3, output, cpus, rank = p
    if read_type == "long":
        a = run(f'{MINIMAP2} ' + params + ' -t '+cpus+' '+refdir+ifile+' '+
                fn1+' '+' > '+fn3 + '  2> ' + output +'logs/mm2log' +
                str(rank) + '.txt',capture_output=True, shell=True)
        assert a.returncode == 0
        
    elif read_type == 'short':
        cmd=f'{BWA} mem '
        if args['simd'] == 'sse': cmd=f'{BWASSE} mem '
        #cmd += f'{BWASSE} mem ' + params + ' -t '+cpus+' '+ refdir+ifile+' '+fn1+' '+' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) + '.txt'
        #else:
        #    cmd=f'{BWA} mem ' + params + ' -t '+cpus+' '+ refdir+ifile+' '+fn1+' '+' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) + '.txt'
        #print('cmd: ', cmd, fn2)
        #print('cmd2: ', refdir, ifile, fn1, fn3, output, rank)
        cmd = cmd + params + ' -t '+cpus+' '+ refdir+ifile+' '+fn1+' ' + fn2 + ' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) + '.txt'
        #print('bwa cmd: ', cmd)
        a=run(cmd, capture_output=True, shell=True)
        assert a.returncode == 0
        
    elif read_type == 'rna_short':
        #cmd = f'{STAR} ' + params + ' --runThreadN '+cpus+' --genomeDir '+ refdir + ' --readFilesIn ' + fn1+' ' + fn2 +' > ' + fn3 + '  2> ' + output +'logs/starlog' + str(rank) + '.txt'
        cmd = f'{STAR} ' + params + ' --runThreadN '+cpus+' --genomeDir '+ refdir + ' --readFilesIn ' + fn1+' ' + fn2 + ' '
        #cmd +=  ' --outFileNamePrefix ' + outfile + str("%05d"%rank) + ' --outTmpDir temp' + str("%05d"%rank)
        cmd += ' --outFileNamePrefix ' + os.path.join(output, 'logs', 'starlog' + str("%05d"%rank))
        cmd += ' --outStd SAM '
        cmd += ' --outTmpDir ' + os.path.join(output, 'temp' + str("%05d"%rank))
        #cmd += '  > '  + output +'logs/starlog' + str(rank) + '.txt'
        cmd += '  > '  + fn3
        #cmd += '  2> ' + output +'logs/starlogerr' + str(rank) + '.txt'
        print('STAR cmd: ', cmd)
        a=run(cmd, capture_output=True, shell=True)
        assert a.returncode == 0
        
    elif read_type == 'meth':
        cmd = "bwameth.py " + params + ' -t ' + cpus + ' --reference ' + refdir+ifile + ' ' + fn1 + ' ' + fn2 + ' > ' + fn3 + ' '
        cmd += ' 2> ' + output +'logs/bwamethlog' + str(rank) + '.txt'
        print('bwa-meth cmd: ', cmd)
        a=run(cmd, capture_output=True, shell=True)

        assert a.returncode == 0

    

#def main(argv):
def main(args):
    global chromo_dict    
    read_type=args["read_type"]        
    ifile=args["refindex"]
    params = args["params"]
    #params= ""
    #if args["params"] != "" and args["params"]  != "None" or args["params"] != None:
    #    params = " -R " +  "\"" + args["params"] + "\""

    read1 = rfile1 = args["read1"]
    read2 = rfile2 = args["read2"]
    read3 = args["read3"]
    readi1 = args["readi1"]
    se_mode = False
    if read2 == "" or read2 == "None":
        se_mode = True
        
    if read_type == "long":
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
    #print('rindex: ', rindex)
    inputdir=args["input"] + "/"
    output=args["output"] + "/"
    tempdir=args["tempdir"]
    if tempdir=="": tempdir=output
    else: tempdir=tempdir + "/"
    refdir=args["refdir"] + "/"
    #if refdir=="": refdir=inputdir
    prof=args["profile"]
    global keep
    #keep=False
    if not args["not_keep_unmapped"]: keep=True
    dindex=args["dindex"]
    
    sample_id=args['sample_id']
    if sample_id == "": sample_id = output
    output_format = args["output_format"]
    
    prefix=args["prefix"]  ## prefix for mutlifq2sortedbam mode reading 'fqprocess' processed fastq files
    suffix=args["suffix"]  ## suffix for mutlifq2sortedbam mode reading 'fqprocess' processed fastq files
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()

    args["rank"], args['nranks'] = rank, nranks

    global ncpus
    ncpus = int(cpus)
    
    if rank == 0:
        if se_mode:
            print("[Info] Running in SE reads only mode w/ ", read_type, ' read type')
        else: print("[Info] Running in PE reads mode!!!")
        if read_type == "short":
            print("[Info] BWA-MEM2 parameters: ", params)
        if read_type == "long":
            print("[Info] mm2-fast parameters: ", params)
    

    def faidx(refdir, ifile):
        g = refdir + ifile
        #print(g.split(".")[-1])
        if g.split(".")[-1] == "gz":
            print("Error: genome file should not be gzip for indexing.\nExiting...")
            return 1
        
        tic=time.time()
        a=run(f'{SAMTOOLS} faidx '+os.path.join(refdir,ifile) + ' > ' + output + 'logs/samlog.txt', capture_output=True, shell=True)
        assert a.returncode==0,"samtools reference index creation failed"
        end=time.time()
        print("\n[Info] Reference to .fai index creation time",end-tic)
        return 0
    #start0 = time.time()
    # Generate data (reads) index file(s)

    ##############################################################################
    flg=0
    if rank == 0 and  dindex == True:
        flg = faidx(refdir, ifile)
        if flg:  
            print("[Info] Exiting due to a failure in creating .fai file.")                    
    allexit(comm, flg)

    ##############################################################################    
    flg = 0
    if rank == 0:
        fo = os.path.join(output, "logs")
        flg = create_folder(fo)
    allexit(comm, flg)        
    ##############################################################################
    
    # Preindex refernce genome if requested
    flg = 0
    if rank==0 and rindex == 'True':
        print("[Info] Indexing Starts", flush=True)
        begin = time.time()
        if read_type == 'meth':
            cmd='bwameth.py index-mem2 ' + refdir+ifile+ ' > ' + output + \
                '/logs/bwamethindexlog.txt'
        elif read_type == 'short':
            cmd = f'{BWA} index '+ refdir + ifile + ' > ' + output + \
                '/logs/bwaindexlog.txt'
        else:
            print("[Error] Only capable of creating index for bwa-mem2 and bwa-meth")
            flg=1
            #sys.exit(1)
        if flg == 0:
            flg = run(cmd, capture_output=True,shell=True)
            flg = flg.returncode
        end=time.time()
        print("\n[Info] bwa index creation time:", end - begin)
        if flg:  print("[Info] bwa index creation failed. Check logs.")
    allexit(comm, flg)                    
    
    ##############################################################################
    input_check(rank, comm, refdir, inputdir, output, se_mode, ifile, rfile1, rfile2, read_type)
    ##############################################################################
        
    # chromo_dict information added
    if keep == False:
        chromo_file = os.path.join(refdir,ifile)+".fai"
        if os.path.isfile(chromo_file):
            f = open(chromo_file, "r")
            lines=f.readlines()
            f.close()
            for line in lines[0:25]:
                chromo=line.split("\t")[0]
                chromo_dict[chromo]= 1
            #if rank == 0:
            #    print(f'chromo_dict: {chromo_dict}')
        else:
            if rank == 0:
                print(f'Error: File "{chromo_file}" not found, resetting keep=True')
            #exit()
            keep = True                


    comm.barrier()                
    #keep = comm.bcast(keep, root=0)        
    #############################################################################
    
    if mode == 'flatmode':
        r'''
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
        
        if rank == 0 and read_type == "short":
            print("[Info] bwa-mem2 starts in flatmode")
        if rank == 0 and read_type == "long":
            print("[Info] mm2-fast starts in flatmode")
            
        # Execute bwamem2 -- may include sort, merge depending on mode
        begin0 = time.time()
        fn1 = pragzip_reader( comm, int(cpus), inputdir+rfile1, output )
        #fn1 = if_reader( comm, int(cpus), inputdir+rfile1, output )
        fn2 = ""
        if not se_mode:
            fn2 = pragzip_reader( comm, int(cpus), inputdir+rfile2, output, last=True )
            #fn2 = if_reader( comm, int(cpus), inputdir+rfile2, output, last=True )
        fn3 = os.path.join(output, outfile + str("%05d"%rank) + ".sam")
        begin1 = time.time()

        p = os.path.join(output,outfile), read_type, params, refdir, ifile, fn1, fn2, fn3, output, cpus, rank
        alignment(p)
        
        end1b=time.time()
        #thr.join()
        comm.barrier()
        end1=time.time()

        if rank==0:
            print("\n[Info] FASTQ to SAM time (flatmode):",end1-begin1)
            print("   (includes wait time:",end1-end1b,")")
            print("[Info] Output SAM files are present in ", output, " folder")

        if rank == 0:
            clean(output)
            
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
            bash_command = '''zcat '''+  inputdir + "/" + barcode_index + \
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
                r1+="--R1 " + inputdir + "/" + rl1[r] + ' '
                r2+="--R2 " + inputdir + "/" + rl2[r] + ' '
                r3+="--R3 " + inputdir + "/" + rl3[r] + ' '

            for r in range(len(ri1)):
                i1+="--I1 " + inputdir + "/" + ri1[r] + ' '

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
                fn1 = os.path.join(inputdir, prefix + "_" + str(r) + ".R1." + suffix)
                fn2 =  os.path.join(inputdir, prefix + "_" + str(r) + ".R3." + suffix)

                if os.path.isfile(fn1) == False or os.path.isfile(fn2) == False:
                    print(f"[Info] Error: Number of files fastq files ({r}) < number of ranks ({nranks})")
                    print(f"!!! Fastq file(s) are not available for processing by rank {r} aborting..\n\n")
                    #sys.exit(1)
                    error_code = 1
                    comm.Abort(error_code)

        comm.barrier()
        #fn1 = "fastq_R1_" + str(rank) + ".fastq.gz"
        #fn2 = "fastq_R3_" + str(rank) + ".fastq.gz"
        fn1 = os.path.join(inputdir, prefix + "_" + str(rank) + ".R1." + suffix)
        fn2 = os.path.join(inputdir, prefix + "_" + str(rank) + ".R3." + suffix)

        #fn1 = inputdir + "/fastq_R1_" + str(rank) + ".fastq.gz"
        #fn2 = inputdir + "/fastq_R3_" + str(rank) + ".fastq.gz"
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
            a=run(f'{BWA} mem ' + params + ' -t '+cpus+' '+refdir+ifile+' '+fn1+' '+
                  fn2+' > '+fn3 + '  2> ' + output +'logs/bwalog' + str(rank) +
                  '.txt',capture_output=True, shell=True)
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

            print("\n[Info] SAM to sort-BAM starts")
            gin2=time.time()

            
    elif mode == 'sortedbam':
        # Preindex refernce genome if requested
        r'''
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
        #print('in read: ', inputdir + rfile1)
        #print("output: ", output)
        fn1 = pragzip_reader( comm, int(cpus), inputdir + rfile1, output )
        #fn1 = if_reader( comm, int(cpus), inputdir + rfile1, output )
        fn2 = ""
        if not se_mode:
            fn2 = pragzip_reader( comm, int(cpus), inputdir + rfile2, output, last=True )
            #fn2 = if_reader( comm, int(cpus), inputdir + rfile2, output, last=True )
        fn3, thr = sam_writer( comm, output+'/aln' )
        #fn3 = output + "/aln" + str(rank) + ".sam"
        begin1 = time.time()
        p = os.path.join(output,outfile), \
            read_type, params, refdir, ifile, fn1, fn2, fn3, output, cpus, rank
        alignment(p)
        end1b = time.time()
        thr.join()
        comm.barrier()
        end1=time.time()

        if rank==0:
            print("\n[Info] FASTQ to SAM time:",end1-begin1)
            print("   (includes wait time:",end1-end1b,")")

            print("\n[Info] SAM to sort-BAM starts")
            begin2 = time.time()

    else:
        print("Error: Incorrect mode specified!")
        print("please specify mode: sortedbam/flatmode/multifq/fqprocessonly")
        os.sys.exit(1)

    # Finish sort, merge, convert to bam depending on mode
    cmd=""
    for i in range(bins_per_rank): 
        binstr = '%05d'%(nranks*i+rank)
        cmd+=f'{SAMTOOLS} sort --threads '+threads+' -T '+tempdir+ '/aln'+binstr+ \
            '.sorted -o '+ output +'/aln'+binstr+'.bam '+ output+'/aln'+binstr+'.sam' + '  2> ' + output +'logs/samsortlog' + str(rank) + '.txt'
        
        if i%20==0:
            a=run(cmd,capture_output=True,shell=True)
            cmd=""
            
    if not cmd=="":
        a=run(cmd,capture_output=True,shell=True)
        assert a.returncode==0,"samtools sort Failed"    
    comm.barrier()

    if rank==0:
        end2=time.time()
        print("[Info] SAM to sort-BAM time:", end2-begin2)

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
        if outfile == "" or outfile == None or outfile == "None":
            outfile = "final_fq2bam"
        #print('outfile: ',outfile)
        cmd+=f'{SAMTOOLS} cat -o ' + os.path.join(output, outfile) + '.sorted.bam ' + infstr
        a=run(cmd,capture_output=True, shell=True)
        assert a.returncode == 0
        print("[Info] Concat done, time taken (s): {:.4f}".format(time.time() - tic))
        #tic2 = time.time()
        os.system('rm ' + output + "/aln*.bam")
        #print("[Info] Concat done, time taken (s): ", time.time() - tic, time.time() - tic2)

        #if rank == 0:
        #    clean_all(output, tempdir)
    elif rank == nranks - 1:
        #print("[Info] cleaning up...")
        tic = time.time()
        if not args["keep_sam"]:
            os.system('rm ' + output + "/*.sam")
        os.system('rm ' + output + "*.idx")
        toc = time.time()
        print("[Info] Clean up done.")
        #print("[Info] cleanup done, time taken (s): ", toc -tic)
        

            
def clean_all(output, tempdir):
    print("[Info] Cleaning up...")
    tic = time.time()
    os.system('rm ' + output + "/*.sam")
    os.system('rm ' + output + "/aln*.bam")
    os.system('rm ' + output + "*.idx")
    toc = time.time()        
    print("[Info] Clean up done, time taken (s): ", toc -tic)

    
def clean(output):
    print("[Info] Cleaning up...")
    tic = time.time()
    os.system('rm ' + output + "*.idx")
    toc = time.time()        
    print("[Info] Clean up done, time taken (s): ", toc -tic)
        
        
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
    #print('args: ', sys.argv)
    #print("finish")
    #main(sys.argv[1:])
    args = json.loads(sys.argv[1])
    #print('argss: ', args)
    main(args)
