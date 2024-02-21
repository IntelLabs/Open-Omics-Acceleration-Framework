#!/usr/bin/python
from subprocess import Popen, PIPE, run
import subprocess
import time
import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def main(argv):
    parser=ArgumentParser()
    parser.add_argument('--input',help="Input data directory")
    parser.add_argument('--output',help="Output data directory")
    parser.add_argument("-i", "--index", help="name of index file")
    parser.add_argument("-r", "--reads", nargs='+',help="name of reads file seperated by space")
    parser.add_argument("-c", "--cpus",default=72,help="Number of cpus. default=72")
    parser.add_argument("-t", "--threads",default=72,help="Number of threads used in samtool operations. default=72")
    #parser.add_argument("-t", "--threads",default=72,help="Number of threads used in samtool operations. default=72")
    parser.add_argument('-in', '--istart',action='store_true',help="It Will start indexing")
    parser.add_argument('-sindex',action='store_true',help="It Will start creating .fai file")
    parser.add_argument('--shards',default=1,help="Number of shards for deepvariant")
    args = vars(parser.parse_args())
    ifile=args["index"]
    rfile1=args["reads"][0]
    rfile2=args["reads"][1]
    cpus=args["cpus"]
    threads=args["threads"]
    index=args["istart"]
    sindex=args["sindex"]
    nproc=args["shards"]
    folder=args["input"]
    output=args["output"]

    t0=time.time()
    file_size = os.path.getsize(folder+rfile1)
    print("\nSize of FASTQ file:",file_size)

    if index==True :
        print("Indexing Starts")
        begin = time.time()
        a=run('../../applications/bwa-0.7.17/bwa index '+folder+ifile,capture_output=True,shell=True)
        end=time.time()
        file_size = os.path.getsize(folder+rfile1)
        print("\nIndex time:",end-begin)
        print("\nSize of FASTQ file:",file_size)
        

    print("bwa starts")
    begin1 = time.time()
    print('../../applications/bwa-0.7.17/bwa mem -t '+cpus+' '+folder+ifile+' '+folder+rfile1+' '+folder+rfile2+' > '+output+'aln.sam')
    a=run('../../applications/bwa-0.7.17/bwa mem -t '+cpus+' '+folder+ifile+' '+folder+rfile1+' '+folder+rfile2+' > '+output+'aln.sam',capture_output=True, shell=True)
    end1=time.time()
    #file_size = os.path.getsize(output+'aln.sam')
    print("\nFASTQ to SAM time:",end1-begin1)
    print("\nSize of SAM file:",file_size)
    
    print("sam to sort-bam starts")
    begin2=time.time()
    print(output+'aln.bam')
    a=run('../../applications/samtools/samtools sort --threads '+threads+' -T /tmp/aln.sorted -o '+output+'aln.bam '+output+'aln.sam',capture_output=True,shell=True)
    end2=time.time()
    file_size = os.path.getsize(output+'aln.bam')
    print("\nSAM to sort-BAM time:",end2-begin2)
    print("\nSize of sort-BAM file",file_size)
    
    begin3=time.time()
    print("Indexing of ref and read starts")
    if sindex==True :
        a=run('../../applications/samtools/samtools faidx '+folder+ifile,capture_output=True,shell=True)
    
    print('../../applications/samtools/samtools index -M -@ '+threads+' '+output+'aln.bam') 
    a=run('../../applications/samtools/samtools index -M -@ '+threads+' '+output+'aln.bam',capture_output=True,shell=True)

    end3=time.time()
    print("\nIndex creation time",end3-begin3)
    
    begin5=time.time()
    #original
    command='sudo docker run -v '+folder+':"/input" -v '+output+':"/output" google/deepvariant:1.5.0 /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/input/'+ifile+' --reads=/output/aln.bam --output_vcf=/output/output.vcf.gz --output_gvcf=/output/output.g.vcf.gz --intermediate_results_dir /output/intermediate_results_dir --num_shards='+nproc+' --dry_run=false'
    #updated
    #command='podman run -v '+folder+':"/input" -v '+output+':"/output" localhost/deepvariant:latest /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/input/'+ifile+' --reads=/output/aln.sorted.new.bam --output_vcf=/output/output.vcf.gz --intermediate_results_dir /output/intermediate_results_dir --num_shards='+nproc+' --pcl_opt --dry_run=false'
    print(command)
    a=run( command+" 2>&1 >> "+output+"log_deepvariant.txt", shell=True)
    #pid=subprocess.call(command,shell=True)
    #os.system(command)
    end5=time.time()
    print("\nDeepVariant runtime",end5-begin5)
    print("Pipeline runtime",end5-t0)

if __name__ == "__main__":
    main(sys.argv[1:])
