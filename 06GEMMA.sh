echo ""
echo ""
echo "Job started on `hostname` at `date`"
echo ""

echo "Prepare soft"
wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
gunzip -c gemma-0.98.5-linux-static-AMD64.gz > gemma-0.98.5-linux-static

echo "Prepare data"
vcftools --gzvcf filter.snps.vcf.gz --max-missing 0.8 --maf 0.05 --max-maf 0.95 --minQ 30 --min-alleles 2 --max-alleles 2 --min-meanDP 10 --max-meanDP 100 --recode --recode-INFO-all --out EpSNP.Miss0.8MAF
vcftools --vcf EpSNP.Miss0.8MAF.vcf --plink --out EpSNP
plink1.9 --file EpSNP --make-bed --out EpSNP
plink1.9 --bfile EpSNP --pca --out EpSNP_covar

echo "run gemma for each pheno"
gemma-0.98.5-linux-static -bfile EpSNP -gk 2 -o EpSNPkin
gemma-0.98.5-linux-static -bfile EpSNP -p EpSNP.pheno.txt -n 1 -k EpSNPkin.sXX.txt -u EpSNP_covar.txt -miss 0.2 -maf 0.05 -lmm 4 -o EpSNP_lmm

echo "Job end on `hostname` at `date`"
