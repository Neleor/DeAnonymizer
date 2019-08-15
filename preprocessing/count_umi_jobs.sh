LOC=/location_you_want_the_split_umi_counts/
OUTLOC=${LOC}/results/
VCF=/path_to_the_genotyped.vcf.gz
BAM=/path_to_cellranger_output/outs/possorted_genome_bam.bam
JOBLOC=${LOC}/jobs/

rm -rf ${JOBLOC}
mkdir ${JOBLOC}
mkdir ${OUTLOC}

chromLengths=(249250621   243199373   198022430   191154276   180915260   171115067   159138663   146364022   141213431   135534747   135006516   133851895   115169878   107349540   102531392   90354753   81195210   78077248   59128983   63025520   48129895   51304566)
stepSize=1250000

for i in {1..22}; do

     end=$((${stepSize}+${chromLengths[i-1]}))
     for ((j=stepSize;j<${end};j+=stepSize)); do
           k="$(($j-$stepSize+1))" 
           echo '#!/bin/bash' >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --job-name=count_umis_'$i.$j >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --output='${JOBLOC}/count_umis_$i.$j.out >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --error='${JOBLOC}/count_umis_$i.$j.err >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --time=11:59:00' >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --cpus-per-task=1' >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --mem-per-cpu=2gb' >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --nodes=1' >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --open-mode=truncate' >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --export=NONE' >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --get-user-env=30L' >> ${JOBLOC}/chr${i}.${j}.sh
           echo '#SBATCH --constraint=tmp02' >> ${JOBLOC}/chr${i}.${j}.sh
           echo ml Python/2.7.12-foss-2015b >> ${JOBLOC}/chr${i}.${j}.sh
           echo 'python /location_of_following_python_file/count_umis_at_vcf_sites.py --input_vcf_file '$VCF' --bamfile '$BAM' --output_file '$OUTLOC/'chr'${i}.${j}'.txt --vcf_region '$i':'$k'-'$j' >' $OUTLOC/${i}.${j}_out.txt '2>' $OUTLOC/${i}.${j}_err.txt  >> ${JOBLOC}/chr${i}.${j}.sh
           sbatch ${JOBLOC}/chr${i}.${j}.sh
    done
done
