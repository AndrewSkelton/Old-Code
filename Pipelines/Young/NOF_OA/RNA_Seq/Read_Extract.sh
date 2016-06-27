
David : 	chr15:89411447-89411560

Cufflinks : chr15:89411408-89411560


samtools view -h -o Flowcell_B_2_Novel_ACAN_exon_David.sam Flowcell_B_2.bam chr15:89411447-89411560
samtools view -h -o Flowcell_B_6_Novel_ACAN_exon_David.sam Flowcell_B_6.bam chr15:89411408-89411560
samtools view -h -o Flowcell_B_2_Novel_ACAN_exon_Cufflinks.sam Flowcell_B_2.bam chr15:89411447-89411560
samtools view -h -o Flowcell_B_6_Novel_ACAN_exon_Cufflinks.sam Flowcell_B_6.bam chr15:89411408-89411560

samtools view -bS Flowcell_B_2_Novel_ACAN_exon_David.sam > Flowcell_B_2_Novel_ACAN_exon_David.bam
samtools view -bS Flowcell_B_6_Novel_ACAN_exon_David.sam > Flowcell_B_6_Novel_ACAN_exon_David.bam
samtools view -bS Flowcell_B_2_Novel_ACAN_exon_Cufflinks.sam > Flowcell_B_2_Novel_ACAN_exon_Cufflinks.bam
samtools view -bS Flowcell_B_6_Novel_ACAN_exon_Cufflinks.sam > Flowcell_B_6_Novel_ACAN_exon_Cufflinks.bam

samtools index Flowcell_B_2_Novel_ACAN_exon_David.bam
samtools index Flowcell_B_6_Novel_ACAN_exon_David.bam
samtools index Flowcell_B_2_Novel_ACAN_exon_Cufflinks.bam
samtools index Flowcell_B_6_Novel_ACAN_exon_Cufflinks.bam

bamToBed -i  Flowcell_B_2_Novel_ACAN_exon_David.bam > Flowcell_B_2_Novel_ACAN_exon_David.bed
bamToBed -i  Flowcell_B_6_Novel_ACAN_exon_David.bam > Flowcell_B_6_Novel_ACAN_exon_David.bed
bamToBed -i  Flowcell_B_2_Novel_ACAN_exon_Cufflinks.bam > Flowcell_B_2_Novel_ACAN_exon_Cufflinks.bed
bamToBed -i  Flowcell_B_6_Novel_ACAN_exon_Cufflinks.bam > Flowcell_B_6_Novel_ACAN_exon_Cufflinks.bed

awk '$2 >= 89411447 && $3 <= 89411560' Flowcell_B_2_Novel_ACAN_exon_David.bed > Flowcell_B_2_Novel_ACAN_exon_David_Filtered.bed
awk '$2 >= 89411447 && $3 <= 89411560' Flowcell_B_6_Novel_ACAN_exon_David.bed > Flowcell_B_6_Novel_ACAN_exon_David_Filtered.bed
awk '$2 >= 89411408 && $3 <= 89411560' Flowcell_B_2_Novel_ACAN_exon_Cufflinks.bed > Flowcell_B_2_Novel_ACAN_exon_Cufflinks_Filtered.bed
awk '$2 >= 89411408 && $3 <= 89411560' Flowcell_B_6_Novel_ACAN_exon_Cufflinks.bed > Flowcell_B_6_Novel_ACAN_exon_Cufflinks_Filtered.bed

bedtools intersect -f 1.0 -b Flowcell_B_2_Novel_ACAN_exon_David_Filtered.bed -abam Flowcell_B_2_Novel_ACAN_exon_David.bam > Flowcell_B_2_Novel_ACAN_exon_David_Intersect.bam
bedtools intersect -f 1.0 -b Flowcell_B_6_Novel_ACAN_exon_David_Filtered.bed -abam Flowcell_B_6_Novel_ACAN_exon_David.bam > Flowcell_B_6_Novel_ACAN_exon_David_Intersect.bam
bedtools intersect -f 1.0 -b Flowcell_B_2_Novel_ACAN_exon_Cufflinks_Filtered.bed -abam Flowcell_B_2_Novel_ACAN_exon_Cufflinks.bam > Flowcell_B_2_Novel_ACAN_exon_Cufflinks_Intersect.bam
bedtools intersect -f 1.0 -b Flowcell_B_6_Novel_ACAN_exon_Cufflinks_Filtered.bed -abam Flowcell_B_6_Novel_ACAN_exon_Cufflinks.bam > Flowcell_B_6_Novel_ACAN_exon_Cufflinks_Intersect.bam

samtools index Flowcell_B_2_Novel_ACAN_exon_David_Intersect.bam
samtools index Flowcell_B_6_Novel_ACAN_exon_David_Intersect.bam
samtools index Flowcell_B_2_Novel_ACAN_exon_Cufflinks_Intersect.bam
samtools index Flowcell_B_6_Novel_ACAN_exon_Cufflinks_Intersect.bam

genomeCoverageBed -split -bg -ibam Flowcell_B_2_Novel_ACAN_exon_David.bam -g hg19.chrom.sizes > Flowcell_B_2_Novel_ACAN_exon_David.bedgraph
genomeCoverageBed -split -bg -ibam Flowcell_B_6_Novel_ACAN_exon_David.bam -g hg19.chrom.sizes > Flowcell_B_6_Novel_ACAN_exon_David.bedgraph
genomeCoverageBed -split -bg -ibam Flowcell_B_2_Novel_ACAN_exon_Cufflinks.bam -g hg19.chrom.sizes > Flowcell_B_2_Novel_ACAN_exon_Cufflinks.bedgraph
genomeCoverageBed -split -bg -ibam Flowcell_B_6_Novel_ACAN_exon_Cufflinks.bam -g hg19.chrom.sizes > Flowcell_B_6_Novel_ACAN_exon_Cufflinks.bedgraph

genomeCoverageBed -split -bg -ibam Flowcell_B_2_Novel_ACAN_exon_David_Intersect.bam -g hg19.chrom.sizes > Flowcell_B_2_Novel_ACAN_exon_David_Intersect.bedgraph
genomeCoverageBed -split -bg -ibam Flowcell_B_6_Novel_ACAN_exon_David_Intersect.bam -g hg19.chrom.sizes > Flowcell_B_6_Novel_ACAN_exon_David_Intersect.bedgraph
genomeCoverageBed -split -bg -ibam Flowcell_B_2_Novel_ACAN_exon_Cufflinks_Intersect.bam -g hg19.chrom.sizes > Flowcell_B_2_Novel_ACAN_exon_Cufflinks_Intersect.bedgraph
genomeCoverageBed -split -bg -ibam Flowcell_B_6_Novel_ACAN_exon_Cufflinks_Intersect.bam -g hg19.chrom.sizes > Flowcell_B_6_Novel_ACAN_exon_Cufflinks_Intersect.bedgraph
