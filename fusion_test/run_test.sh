#!/usr/bin/env bash

outdir=tophat_out_fusion

for file in test_*.fasta
do
   rm -rf $outdir ; ~/work/tophat/tophat_trunk/src/tophat -o $outdir --fusion-do-not-resolve-conflicts --max-intron-length 500 --fusion-min-dist 500 --keep-tmp --fusion-search --bowtie1 --keep-fasta-order testcases/test $file >& /dev/null
   echo "$file:"
   total=`grep ">" $file | wc -l`
   found=`samtools view $outdir/accepted_hits.bam | awk '{print $1}' | sort | uniq | wc -l`
   echo "$found/$total"
done
