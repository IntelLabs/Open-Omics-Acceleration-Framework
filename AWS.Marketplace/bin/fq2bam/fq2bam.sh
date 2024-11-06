docker run -v ~/data/fq2bam/input:/input -v  ~/data/fq2bam/output:/out -v  ~/data/fq2bam/refdir:/refdir -v  ~/data/fq2bam/tempdir:/tempdir fq2bam:latest bash run_bwa.sh sortedbam /input/config.yaml

