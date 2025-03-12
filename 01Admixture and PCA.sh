#filter
vcftools --gzvcf filter.snps.vcf.gz --max-maf 0 --minQ 30 --max-missing 0.8 --min-meanDP 10 --max-meanDP 100 --recode --stdout | bgzip -c > invariant.vcf.gz 
vcftools --gzvcf filter.snps.vcf.gz --mac 1 --min-alleles 2 --max-alleles 2 --minQ 30 --max-missing 0.8 --min-meanDP 10 --max-meanDP 100 --recode --stdout | bgzip -c > variant.vcf.gz  
tabix invariant.vcf.gz
tabix variant.vcf.gz
bcftools concat --allow-overlaps variant.vcf.gz invariant.vcf.gz -O z -o allsites.vcf.gz

#LD 
/home/ZhuXL/Biosoft/plink1.9/plink --vcf variant.vcf.gz --indep-pairwise 50 5 0.2 \
--out CPA_LD --allow-extra-chrzcat QMsnp_MS90_MAF.vcf.gz | awk -v target_file="CPA_LD.prune.in" \
-F '\t' 'BEGIN {while(getline < target_file) targets[$1]=1} /^#/ || $3 in targets {print}' | bgzip -c > CPA_LD.vcf.gz 
#Admixture
plink --allow-extra-chr --vcf CPA_LD.vcf.gz --make-bed --out CPA.ADM
#admixture
for K in {1..20}; do admixture --cv CPA.ADM $K -j20 |tee log${K}.out;done 

#PCA
VCF2PCACluster -InVCF CPA_LD.vcf.gz -KinshipMethod 4 -OutPut PCA 