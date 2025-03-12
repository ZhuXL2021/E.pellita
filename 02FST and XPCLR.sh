#Pi
vcftools --gzvcf CPA_LD.vcf.gz  --keep popmap_Pop1.txt --window-pi 10000  --out Pop1
vcftools --gzvcf CPA_LD.vcf.gz  --keep popmap_Pop2.txt --window-pi 10000  --out Pop2
vcftools --gzvcf CPA_LD.vcf.gz  --keep popmap_Pop3.txt --window-pi 10000  --out Pop3
#FST；
vcftools --gzvcf CPA_LD.vcf.gz --weir-fst-pop popmap_Pop1.txt --weir-fst-pop popmap_Pop2.txt --fst-window-size 10000 --fst-window-step  10000 --out Pop1_Pop2  
vcftools --gzvcf CPA_LD.vcf.gz --weir-fst-pop popmap_Pop1.txt --weir-fst-pop popmap_Pop3.txt --fst-window-size 10000 --fst-window-step  10000 --out Pop1_Pop3
vcftools --gzvcf CPA_LD.vcf.gz --weir-fst-pop popmap_Pop2.txt --weir-fst-pop popmap_Pop3.txt --fst-window-size 10000 --fst-window-step  10000 --out Pop2_Pop3

#LDdecay
vcftools --gzvcf CPA_LD.vcf.gz --recode --recode-INFO-all --stdout --keep popmap_Pop1.txt | bgzip -c > 15Pop1.vcf.gz 
vcftools --gzvcf CPA_LD.vcf.gz --recode --recode-INFO-all --stdout --keep popmap_Pop2.txt | bgzip -c > 15Pop2.vcf.gz 
vcftools --gzvcf CPA_LD.vcf.gz --recode --recode-INFO-all --stdout --keep popmap_Pop3.txt | bgzip -c > 15Pop3.vcf.gz 
PopLDdecay -InVCF 15Pop1.vcf.gz -OutStat 15Pop1
PopLDdecay -InVCF 15Pop2.vcf.gz -OutStat 15Pop2
PopLDdecay -InVCF 15Pop3.vcf.gz -OutStat 15Pop3


#XPCLR
for i in {1..11}; do
  echo $i
done | xargs -n 1 -P 11 -I {} sh -c '
xpclr -I CPA_LD.vcf.gz -O ./Pop1_Pop3.Chr{}.xpclr-python.out -C {} -Sa ./Pop1.list -Sb ./Pop3.list --ld 0.95 --maxsnps 200 --size 10000 \
    --step 10000 --rrate 1e-8
'