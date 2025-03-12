echo ""
echo ""
echo "Job started on `hostname` at `date`"
echo ""

echo "Prepare soft"
wget https://github.com/BioH4z/gatk3.8.0/blob/master/GenomeAnalysisTK.jar
#Install Java
sudo mkdir -p /usr/local/java
cd /usr/local/java
sudo wget -c --header http://download.oracle.com/otn-pub/java/jdk/8u131-b11/d54c1d3a095b4ff2b6607d096fa80163/jdk-8u131-linux-x64.tar.gz
sudo tar xvzf jdk-8u131-linux-x64.tar.gz
echo 'export PATH=$PATH:/usr/local/java/jdk1.8.0_131/jre/bin' >> ~/.bashrc
source ~/.bashrc
#Download Picard tools
wget https://github.com/broadinstitute/picard/releases/download/3.0.0/picard.jar

echo "prepare genome reference"
bwa index -a bwtsw Egrandis_297_v2.0.fa
java -Xmx20g -jar picard.jar CreateSequenceDictionary R=Egrandis_297_v2.0.fa O=Egrandis_297_v2.0.dict
samtools faidx Egrandis_297_v2.0.fa

#The faster way to analysis each sample using screen, not using the below command
echo "Add Read group information and do mapping"
 for R1 in *_R1.fq.gz;do
    SM=$(echo $R1 | cut -d"_" -f1)                                          ##sample ID
    LB=$(echo $R1 | cut -d"_" -f1,2)                                        ##library ID
    PL="Illumina"                                                           ##platform (illumina)
    RGID=$(zcat $R1 | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)       ##read group identifier 
    PU=$RGID.$LB                                                            ##Platform Unit
    echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"

    R2=$(echo $R1 | sed 's/_R1_/_R2_/')
    echo $R1 $R2
    bwa mem -t 4 -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" Egrandis_297_v2.0.fa $R1 $R2 > ${R1%_R1_001.pe.fq.gz}.sam
  done
  
echo "Sort BAM file by coordinate, convert to BAM"
for samfile in *.sam;do
  sample=${samfile%.sam}
  samtools view -bS -o $sample.bam $samfile
  samtools sort $sample.bam $sample.sorted
done

echo "Mark Duplicates"
for sample in *.sorted.bam;do
  name=${sample%.sorted.bam}
  java -Xmx30g -jar picard.jar MarkDuplicates INPUT=$sample OUTPUT=$name.dedup.bam METRICS_FILE=$name.metrics.txt;
done

echo "Build BAM Index"
for sample in *.dedup.bam;do
  name=${sample%.dedup.bam}
  java -Xmx30g -jar picard.jar BuildBamIndex INPUT=$sample
done

echo "Call Variants"
for sample in *.dedup.bam;do
  name=${sample%.dedup.bam}
java -Xmx30g -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R Egrandis_297_v2.0.fa -I $sample --emitRefConfidence GVCF -o $name.raw.vcf -nct 2
done

echo "Joint Genotyping"
ls *.vcf > input_gvcfs.list
java -Xmx80g -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R Egrandis_297_v2.0.fa -V input_gvcfs.list -stand_call_conf 10 -ploidy 2 -maxAltAlleles 2 -o raw.GVCFall.vcf

echo "Filters Genotyping"
java -Xmx30g -jar GenomeAnalysisTK.jar -T SelectVariants -R Egrandis_297_v2.0.fa -V raw.GVCFall.vcf -selectType SNP -o raw.SNPs.vcf
java -Xmx30g -jar GenomeAnalysisTK.jar -T SelectVariants -R Egrandis_297_v2.0.fa -V raw.GVCFall.vcf -selectType INDEL -o raw.INDELs.vcf

gatk --java-options "-Xmx80g" VariantFiltration -R Egrandis_297_v2.0.fa -V raw.SNPs.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "SNP_FILTER" -O filter.snps.vcf.gz

gatk --java-options "-Xmx80g" VariantFiltration -R Egrandis_297_v2.0.fa -V raw.INDELs.vcf --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -20.0" --filter-name "INDEL_FILTER" -O filter.indels.vcf.gz

echo "Job end on `hostname` at `date`"

