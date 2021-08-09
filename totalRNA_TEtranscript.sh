#!/bin/bash

module load gcc/6.2.0
module load python/2.7.12
module load samtools/1.9
module load stringtie/1.3.3b

source ~/.bash_profile

GTF=/n/groups/zhanglab/zychen/genomes_annotations/mm10/annotations/mm10_plus_ercc92_fixed.gtf
REPEAT=/n/groups/zhanglab/zychen/genomes_annotations/mm10/annotations/mm10_rmsk_TE.gtf

TEtranscripts -t sample1.hisat2.bam sample2.hisat2.bam \
              -c sample3.hisat2.bam sample4.hisat2.bam \
                 --GTF $GTF --TE $REPEAT --stranded reverse \
	         >& TEtranscript.log

exit;

