#!/bin/bash

samtools merge -@4 --write-index WT.cleanAligned.sortedByCoord.out.markDups.bam A171_WT1_1.cleanAligned.sortedByCoord.out.markDups.bam A171_WT2_1.cleanAligned.sortedByCoord.out.markDups.bam A171_WT3_1.cleanAligned.sortedByCoord.out.markDups.bam
samtools merge -@4 --write-index M4.cleanAligned.sortedByCoord.out.markDups.bam A1714_1_1.cleanAligned.sortedByCoord.out.markDups.bam A1714_2_1.cleanAligned.sortedByCoord.out.markDups.bam A1714_3_1.cleanAligned.sortedByCoord.out.markDups.bam
samtools merge -@4 --write-index M6.cleanAligned.sortedByCoord.out.markDups.bam A1716_1_1.cleanAligned.sortedByCoord.out.markDups.bam A1716_2_1.cleanAligned.sortedByCoord.out.markDups.bam A1716_3_1.cleanAligned.sortedByCoord.out.markDups.bam

bamCoverage --bam WT.cleanAligned.sortedByCoord.out.markDups.bam --normalizeUsing BPM -bs 1 -p 16 -o WT.bw
bamCoverage --bam M4.cleanAligned.sortedByCoord.out.markDups.bam --normalizeUsing BPM -bs 1 -p 16 -o M4.bw
bamCoverage --bam M6.cleanAligned.sortedByCoord.out.markDups.bam --normalizeUsing BPM -bs 1 -p 16 -o M6.bw

