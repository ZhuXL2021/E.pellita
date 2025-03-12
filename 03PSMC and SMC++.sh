#PSMC
for i in `cat psmc.list`; do
  echo $i
done | xargs -n 1 -P 15 -I {} sh -c '
samtools mpileup -C50 -uf ref.fasta {}.sorted.dd.header.bam | \
bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 200| gzip > {}.psmc.fq.gz
fq2psmcfa -q20 {}.psmc.fq.gz > {}.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o {}.psmc {}.psmcfa
perl psmc_plot.pl -u 4.77e-9 -g 25 -p {}_psmc.plot {}.psmc
' 
#bootstrap
cat psmc.list | while read i; do
    seq 10 | xargs -I {} sh -c 'psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o ./round10/${1}.round-{}.psmc ${1}.psmcfa' _ "$i" &
done
wait


#plot
perl psmc_plot.pl -u 4.9e-8 -g 5 -p -X 2000000 -Y 5 mean.plot all.psmc

#SMC++
bwa index genome.fa
samtools faidx genome.fa
picard -Xmx4G CreateSequenceDictionary R=genome.fa O=genome.dict
mkdir -p Pop1 Pop2 Pop3
for i in {1..11}
do
smc++ vcf2smc --mask genome.mask.bed.gz CPA_LD.vcf.gz Pop1/Chr${i}.smc.gz ${i} Pop1:S-118,S-121,S-123,S-125,S-218,S-245,S-248,S-253,S-26,S-37,S-38,S-40,S-42,S-43,S-49
smc++ vcf2smc --mask genome.mask.bed.gz CPA_LD.vcf.gz Pop2/Chr${i}.smc.gz ${i} Pop2:S-137,S-159,S-190,S-195,S-256,S-258,S-259,S-260,S-262,S-264,S-271,S-272,S-275,S-276,S-279
smc++ vcf2smc --mask genome.mask.bed.gz CPA_LD.vcf.gz Pop3/Chr${i}.smc.gz ${i} Pop3:S-103,S-108,S-230,S-232,S-235,S-238,S-239,S-240,S-241,S-242,S-243,S-50,S-93,S-94,S-95
done
smc++ estimate --spline cubic --knots 30 --cores 10 -o Pop1/ 4.9e-8 Pop1/*.smc.gz 
smc++ estimate --spline cubic --knots 30 --cores 10 -o Pop2/ 4.9e-8 Pop2/*.smc.gz 
smc++ estimate --spline cubic --knots 30 --cores 10 -o Pop3/ 4.9e-8 Pop3/*.smc.gz 
#plot
smc++ plot plot.Pop1.pdf Pop1/*.final.json -x 100 100000 -g 5  -y 0 60000 -c 
smc++ plot plot.Pop2.pdf Pop2/*.final.json -x 100 100000 -g 5  -y 0 60000 -c 
smc++ plot plot.Pop3.pdf Pop3/*.final.json -x 100 100000 -g 5  -y 0 60000 -c 
#split time
mkdir Pop2.Pop3 Pop2.Pop3.split && cd Pop2.Pop3
ls ../Pop2/*.smc.gz|while read line;do ln -s ${line} Pop2$(basename ${line}); done
ls ../Pop3/*.smc.gz|while read line;do ln -s ${line} Pop3$(basename ${line}); done
cd ../
nohup smc++ split -o Pop2.Pop3.split/ Pop2/model.final.json \
Pop3/model.final.json Pop2.Pop3/*.smc.gz --cores 15 > Pop2.Pop3.split.log 2>&1& 
smc++ plot Pop2.Pop3.pdf Pop2.Pop3.split/model.final.json -x 100 100000 -g 5  -y 0 60000 -c


