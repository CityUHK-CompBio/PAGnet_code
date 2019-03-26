#used MACS software (version 2.1.0) to call narrow peaks#
#here is peak calling code for AlgR#
macs2 callpeak -t algR.fastq.bam -c Input.fastq.bam -f BAM -g 6264404 --nomodel --extsize 75 -p 1e-5 -n algR
