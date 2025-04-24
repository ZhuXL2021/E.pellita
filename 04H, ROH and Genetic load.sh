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

#genetic load
#building database
java -jar -Xmx24g snpEff.jar build -c snpEff.config -gff3 -v CPA -d -noCheckCds -noCheckProtein
#ancestral states
vcftools --gzvcf E.urophylla.pass.snps.vcf.gz --min-alleles 2 --max-alleles 2 --minQ 30 --min-meanDP 10 --max-meanDP 200 --recode --stdout | bgzip -c > OGind.vcf.gz
mkdir -p RefAnc AltAnc OG 
vcftools --gzvcf ./OGind.vcf.gz --counts --out ./OG/outgroup
cat ./OG/outgroup.frq.count | awk -F"[:\t]" '{if($4 == $6 && $4 != 0 ) print $1"\t"$2 }' >> ./RefAnc/outgroup.homo.Ref.anc.pos
cat ./OG/outgroup.frq.count | awk -F"[:\t]" '{if($4 == $8 && $4 != 0 ) print $1"\t"$2 }' >> ./AltAnc/outgroup.homo.Alt.anc.pos 
#As Ref.anc for example
mkdir -p ./Cout ./Pos ./IndVCF ./temp ./Output ./stat ./lof
for n in `cat ind.list`; do
  echo $n
done | xargs -n 1 -P 20 -I {} sh -c '
   vcftools --gzvcf CPA_LD.vcf.gz --indv {} --positions outgroup.homo.Ref.Anc.pos --max-missing 1.0 --counts --out "./Cout/{}.Ref"
   '
for n in `cat ind.list`; do
   cat ./Cout/${n}.Ref.frq.count | awk -F "[\t/:]" '{if($8 == 2) print $1"\t"$2}' >>  ./Pos/${n}.Ref.homo.pos ###When the ancestor state is Alt, change $8 ==2 to $6==2
   cat ./Cout/${n}.Ref.frq.count | awk -F "[\t/:]" '{if($8 == 1) print $1"\t"$2}' >>  ./Pos/${n}.Ref.hete.pos ###When the ancestor state is Alt, change $8 ==2 to $6==2 
done
for n in `cat ind.list`; do
  echo $n
done | xargs -n 1 -P 20 -I {} sh -c '  
   vcftools --gzvcf CPA_LD.vcf.gz --indv {} --positions ./Pos/{}.Ref.homo.pos --recode --out ./IndVCF/{}.Ref.homo 
   vcftools --gzvcf CPA_LD.vcf.gz --indv {} --positions ./Pos/{}.Ref.hete.pos --recode --out ./IndVCF/{}.Ref.hete
   '
for n in `cat ind.list`; do   
####annotation
   java -jar -Xmx50g snpEff.jar -c snpEff.config -no-downstream -no-upstream -v CPA ./IndVCF/${n}.Ref.homo.recode.vcf >> ./IndVCF/${n}.Ref.homo.anc
   java -jar -Xmx50g snpEff.jar -c snpEff.config -no-downstream -no-upstream -v CPA ./IndVCF/${n}.Ref.hete.recode.vcf >> ./IndVCF/${n}.Ref.hete.anc
   grep  missense ./IndVCF/${n}.Ref.homo.anc >> ./temp/${n}.Ref.homo.missense.vcf
   grep  missense ./IndVCF/${n}.Ref.hete.anc >> ./temp/${n}.Ref.hete.missense.vcf
###input
   cat ./temp/${n}.Ref.homo.missense.vcf | grep -v '#' | awk -F"\t|\\\\|" '{printf ("%s\t%s\t%s\t%s\t%s\n", $18,$1,$2,$4,$5)}' |  awk -F"\t" '{OFS = FS} { gsub(/p\./,"", $1); print }' | awk '{sub(/[0-9]+/," ",$1); print}'  >> ./temp/${n}.Ref.homo.missense.input
   cat ./temp/${n}.Ref.hete.missense.vcf | grep -v '#' | awk -F"\t|\\\\|" '{printf ("%s\t%s\t%s\t%s\t%s\n", $18,$1,$2,$4,$5)}' |  awk -F"\t" '{OFS = FS} { gsub(/p\./,"", $1); print }' | awk '{sub(/[0-9]+/," ",$1); print}'  >> ./temp/${n}.Ref.hete.missense.input
   done
###DEL  
for n in `cat ind.list`; do   
   python yuli.corrected.py ./temp/${n}.Ref.homo.missense.input  ./temp/${n}.Ref.homo.missense.output.temp
   python yuli.corrected.py ./temp/${n}.Ref.hete.missense.input  ./temp/${n}.Ref.hete.missense.output.temp
   sed -e "s/\s\+/\n/g" ./temp/${n}.Ref.homo.missense.output.temp | awk '{if(NR%3!=0)ORS=" ";else ORS="\n"}1' >> ./Output/${n}.Ref.homo.missense.output
   sed -e "s/\s\+/\n/g" ./temp/${n}.Ref.hete.missense.output.temp | awk '{if(NR%3!=0)ORS=" ";else ORS="\n"}1' >> ./Output/${n}.Ref.hete.missense.output
   done
for n in `cat ind.list`; do
   paste ./temp/${n}.Ref.hete.missense.input ./Output/${n}.Ref.hete.missense.output >./stat/${n}.Ref.hete.missense
   paste ./temp/${n}.Ref.homo.missense.input ./Output/${n}.Ref.homo.missense.output >./stat/${n}.Ref.homo.missense
   done  
