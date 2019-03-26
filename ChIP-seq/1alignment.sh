# !/bin/sh
for id in ~/PA/ChIP-seq/fastq/*fastq
do
echo $id
bowtie -t -k 1 -m 1 -p 20 ~/PA/ChIP-seq/PA_Index $id --sam | samtools view -q 255 -bS - > $id.bam
done
