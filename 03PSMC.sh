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
