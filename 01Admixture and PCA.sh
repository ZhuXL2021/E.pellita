#Admixture
plink --allow-extra-chr --vcf 1286304SNP.vcf.gz --make-bed --out CPA.ADM
#admixture
for K in {1..20}; do admixture --cv CPA.ADM $K -j20 |tee log${K}.out;done 

#PCA
VCF2PCACluster -InVCF 1286304SNP.vcf.gz -KinshipMethod 4 -OutPut PCA 