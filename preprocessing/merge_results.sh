NAME=/location_you_want_the_merged_file/chrALL.molecule_counts.txt
LOC=/location_of_your_results/

head -1 $LOC/chr1.12500000.txt > $NAME

chromLengths=(249250621   243199373   198022430   191154276   180915260   171115067   159138663   146364022   141213431   135534747   135006516   133851895   115169878   107349540   102531392   90354753   81195210   78077248   59128983   63025520   48129895   51304566 )
stepSize=1250000

for i in {1..22}; do
    end=$((${stepSize} + ${chromLengths[i-1]}))

    for ((j=stepSize;j<${end};j+=stepSize)); do
     tail -n +2 $LOC/chr${i}.${j}.txt >> $NAME
done

done
