#!/usr/bin/bash

CALLER=`basename $0`         # The Caller name
usage()
{
        echo "Name: a01mutationCall.sh"
        echo "DESCRIPTION: To detect structure variation by multi-tools in tumor-normal paired samples."
        echo "USAGE: $CALLER [-h] [-s] <FI: ref.fa>  <FI: tumor.bam> <FI:normal.bam>"
        echo "WHERE: -h = help       "
        echo "       -s = silent (no prompts)"
        echo "PREREQUISITES:"
        echo "* The environment variable windowSize ,ref ,map and bam  must be set"
        echo "NOTES: sh ./a01mutationCall.sh ref.fa tumor.bam normal.bam"
        echo "AUTHOR: leimengyue, leimengyue@genomics.cn"
        echo "ORGANIZATION: BGI shenzhen"
        echo "VERSION: 1.0"
        echo "CREATED: 07/10/2017 10:58:27 AM"
        echo "REVISION: ---"
        echo "$CALLER: exiting now ."
        exit 1
}
if [ $# != 3 ]
then
        echo "Environment variable is not set right."
        usage
fi
echo "$CALLER start now " 
echo `date`


ref = $1
tumorBam = $2
normalBam = $3


***************************************************************
#Crest
***************************************************************

perl CREST/extractSClip.pl -i $tumorBam --ref_genome $ref -o Crest/ 
perl CREST/CREST.pl --blatserver cngb-oxcompute-1 --blatport 9002 -f Crest/$tumorBam\.cover -d $tumorBam -g $normalBam --ref_genome $ref -t /hwfssz1/ST_CANCER/POL/SHARE/DataBase/hg19/GATKv2.8_bundle/ucsc.hg19.2bit -o Crest/

***************************************************************
#Delly(DEL, DUP, INV, TRA, INS)
***************************************************************

delly_parallel_linux_x86_64bit call -t INV -g $ref -o Delly/INV.bcf $tumorBam $normalBam
delly_parallel_linux_x86_64bit filter -t INV -f somatic -o Delly/INV.filter.bcf -s samples.tsv -g $ref Delly/INV.bcf
bcftools view Delly/INV.filter.bcf |grep -v \# >Delly/INV.filter.vcf

delly_parallel_linux_x86_64bit call -t INS -g $ref -o Delly/INS.bcf $tumorBam $normalBam
delly_parallel_linux_x86_64bit filter -t INS -f somatic -o Delly/INS.filter.bcf -s samples.tsv -g $ref Delly/INS.bcf
bcftools view Delly/INS.filter.bcf |grep -v \# >Delly/INS.filter.vcf

delly_parallel_linux_x86_64bit call -t TRA -g $ref -o Delly/TRA.bcf $tumorBam $normalBam
delly_parallel_linux_x86_64bit filter -t TRA -f somatic -o Delly/TRA.filter.bcf -s samples.tsv -g $ref Delly/TRA.bcf
bcftools view Delly/TRA.filter.bcf |grep -v \# >Delly/TRA.filter.vcf

delly_parallel_linux_x86_64bit call -t DEL -g $ref -o Delly/DEL.bcf $tumorBam $normalBam
delly_parallel_linux_x86_64bit filter -t DEL -f somatic -o Delly/DEL.filter.bcf -s samples.tsv -g $ref Delly/DEL.bcf
bcftools view Delly/DEL.filter.bcf |grep -v \# >Delly/DEL.filter.vcf

delly_parallel_linux_x86_64bit call -t DUP -g $ref -o Delly/DUP.bcf $tumorBam $normalBam
delly_parallel_linux_x86_64bit filter -t DUP -f somatic -o Delly/DUP.filter.bcf -s samples.tsv -g $ref Delly/DUP.bcf
bcftools view Delly/DUP.filter.bcf |grep -v \# >Delly/DUP.filter.vcf

cat Delly/INV.filter.vcf Delly/INS.filter.vcf Delly/TRA.filter.vcf Delly/DEL.filter.vcf Delly/DUP.filter.vcf > Delly/somatic.vcf

***************************************************************
#Lumpy
***************************************************************

# Extract the discordant paired-end alignments
samtools view -b -F 1294 $normalBam >$normalBam\.discordants.unsorted.bam
# Extract the split-read alignments
samtools view -h $normalBam | lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > $normalBam\.splitters.unsorted.bam
samtools sort $normalBam\.discordants.unsorted.bam $normalBam\.discordants.bam
samtools sort $normalBam\.splitters.unsorted.bam $normalBam\.splitters.bam

# Extract the discordant paired-end alignments
samtools view -b -F 1294 $tumorBam >$tumorBam\.discordants.unsorted.bam
# Extract the split-read alignments
samtools view -h $tumorBam | lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > $tumorBam\.splitters.unsorted.bam
samtools sort $tumorBam\.discordants.unsorted.bam $tumorBam\.discordants.bam
samtools sort $tumorBam\.splitters.unsorted.bam $tumorBam\.splitters.bam

lumpy-sv/bin/lumpyexpress -B $tumorBam,$normalBam -S $tumorBam\.splitters.bam,$normalBam\.splitters.bam -D $tumorBam\.discordants.bam,$normalBam\.discordants.bam -o tumor_normal.vcf

***************************************************************
#Manta
***************************************************************

manta.centos5_x86_64/bin/configManta.py --normalBam $normalBam --tumorBam $tumorBam --referenceFasta $ref --runDir Manta/

***************************************************************
#novoBreak

sh novoBreak_distribution/run_novoBreak.sh novoBreak_distribution/ $ref $tumorBam $normalBam novoBreak/

***************************************************************
