#H
for bam in `cat ind.list`; do

echo -e "angsd: $(date)"
echo "Run angsd on $bam"
outname=ReHo_Mor_SITE_$bam
# -i 
angsd -i ${bam}.dedup.bam -GL 2 -dosaf 1 -anc $ref -ref $ref -doCounts 1 \
     -minQ 20 -minmapq 30 
     -nthreads 20 -out ${outname}
echo -e "theta: $(date)"
echo "Run theta on $bam"
# for global heterozigosity
realSFS ${outname}.saf.idx -P 4 -fold 1 -bootstrap 10 > ${outname}.b10.est.ml
# tole = tolerence_for_breaking_EM
cat ${outname}.b10.est.ml | R -e "f <- file('stdin');a<-scan(f);print(a[2]/sum(a))"
done

#ROH
plink --allow-extra-chr --vcf CPA_LD.vcf.gz --homozyg --homozyg-window-het 0 --homozyg-snp 50 --homozyg-kb 10 \
--homozyg-density 5000 --homozyg-gap 5000 --out ROHoutput

