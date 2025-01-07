#Pi
vcftools --gzvcf 1286304SNP.vcf.gz  --keep ind.txt --window-pi 10000  --out ind.pi

#ROH
plink --allow-extra-chr --vcf 1286304SNP.vcf.gz --homozyg --homozyg-window-het 0 --homozyg-snp 50 --homozyg-kb 10 --homozyg-density 5000 --homozyg-gap 5000 --out ROHoutput

#genetic load
#building database
java -jar -Xmx24g snpEff.jar build -c snpEff.config -gff3 -v CPA -d -noCheckCds -noCheckProtein
#ancestral states
vcftools --gzvcf 1286304SNP.vcf.gz  --hardy  --out ./Pop/population 
awk -F'\t' 'NR>1 && $3 ~ /^([0-9]+\/)+([0-9]+)$/ {split($3, arr, "/"); if(arr[1]>148) print $1"\t"$2}' ./Pop/population.hwe > ./RefAnc/population.homo.Ref.Anc.pos 
awk -F'\t' 'NR>1 && $3 ~ /^[0-9]+\/[0-9]+\/([0-9]+)/ {split($3, arr, "/"); if(arr[3]>148) print $1"\t"$2}'  ./Pop/population.hwe > ./AltAnc/population.homo.Alt.Anc.pos 
#As Ref.anc for example
mkdir -p ./Cout ./Pos ./IndVCF ./temp ./Output ./stat ./lof
for n in `cat ind.list`; do
  echo $n
done | xargs -n 1 -P 20 -I {} sh -c '
   vcftools --gzvcf 1286304SNP.vcf.gz --indv {} --positions population.homo.Ref.Anc.pos --max-missing 1.0 --counts --out "./Cout/{}.Ref"
   '
for n in `cat ind.list`; do
   cat ./Cout/${n}.Ref.frq.count | awk -F "[\t/:]" '{if($8 == 2) print $1"\t"$2}' >>  ./Pos/${n}.Ref.homo.pos ###When the ancestor state is Alt, change $8 ==2 to $6==2
   cat ./Cout/${n}.Ref.frq.count | awk -F "[\t/:]" '{if($8 == 1) print $1"\t"$2}' >>  ./Pos/${n}.Ref.hete.pos ###When the ancestor state is Alt, change $8 ==2 to $6==2 
done
for n in `cat ind.list`; do
  echo $n
done | xargs -n 1 -P 20 -I {} sh -c '  
   vcftools --gzvcf 1286304SNP.vcf.gz --indv {} --positions ./Pos/{}.Ref.homo.pos --recode --out ./IndVCF/{}.Ref.homo 
   vcftools --gzvcf 1286304SNP.vcf.gz --indv {} --positions ./Pos/{}.Ref.hete.pos --recode --out ./IndVCF/{}.Ref.hete
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
###LOF
   cat ./IndVCF/${n}.Ref.homo.anc | java -jar -Xmx3g SnpSift.jar filter "LOF[*].PERC=1.00" >> ./lof/${n}.Ref.homo.LOF.vcf
   cat ./IndVCF/${n}.Ref.hete.anc | java -jar -Xmx3g SnpSift.jar filter "LOF[*].PERC=1.00" >> ./lof/${n}.Ref.hete.LOF.vcf
   awk '$1 !~ /^#/ { print $1, $2, $4, $5 }' ./lof/${n}.Ref.homo.LOF.vcf > ./lof/${n}.Ref.homo.LOF
   awk '$1 !~ /^#/ { print $1, $2, $4, $5 }' ./lof/${n}.Ref.hete.LOF.vcf > ./lof/${n}.Ref.hete.LOF
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
