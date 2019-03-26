# !/bin/sh
for id in ~/PA/alignment/*bam
do
echo $id
samtools sort -o ${id}.sorted.bam $id
samtools index ${id}.sorted.bam
bamCoverage --bam ${id}.sorted.bam -o ${id}.bw --binSize 20 --effectiveGenomeSize 6264404
done
