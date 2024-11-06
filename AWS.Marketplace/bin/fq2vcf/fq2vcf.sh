echo "Running part1/2"
docker run -v /home/ubuntu/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/config:/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/scripts/aws/config -v ~/data/fq2vcf/input:/reads -v ~/data/fq2vcf/refdir:/ref -v ~/data/fq2vcf/output:/output -it deepvariant:part1 bash run_pipeline_ec2_part1.sh

echo "Running part2/2"
docker run -v /home/ubuntu/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/extra_scripts/config:/opt/deepvariant/config  -v ~/data/fq2vcf/refdir:/ref -v ~/data/fq2vcf/output:/output -it deepvariant:part2 bash run_pipeline_ec2_part2.sh

echo "Cleaning"
cd ~/data/fq2vcf/output/
sudo rm -rf ~/data/fq2vcf/output/0000*
sudo rm -rf ~/data/fq2vcf/output/intermediate*
sudo rm bin_region.pkl *.bam *.bai *.idx
cd -
echo "Cleaning done"
echo "Execution completed."
