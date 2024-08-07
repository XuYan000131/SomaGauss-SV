#######################
#minimap2
#! /bin/bash

clean_fq=$1
sampleid=$2
workdir=$3
ref=$4
minimap2=$5
samtools=$6
#inputbam=$5

#clean_fq=$workdir/clean.$sampleid.fq.gz
#clean_fq=$workdir/clean_Q7_L1000_$sampleid.fq.gz
sam=$workdir/${sampleid}.sam
bam=$workdir/${sampleid}.bam
clean_fastq_plots=$workdir/${sampleid}.clean-fastq-plots
bam_sort=$workdir/${sampleid}_minimap2.2.26_sorted.bam
bam_sort_tag=$workdir/${sampleid}_minimap2.2.26_sorted_tag.bam

#echo "$(date) 1. Start QC: $sampleid"
#gunzip -c $fq|/home/wg_fzt/miniconda3/envs/nanoplot/bin/NanoFilt -q 7 -l 1000|gzip > $clean_fq
#echo "$(date) 1. Finish QC: $sampleid"

#echo "$(date) 2. Start NanoPlot: $sampleid"
#/home/wg_fzt/miniconda3/envs/nanoplot/bin/NanoPlot  --fastq $clean_fq \
#   -o $clean_fastq_plots \
#   --maxlength 40000 \
#   -t 40 \
#   --plots hex dot
#echo "$(date) 2. Finish NanoPlot: $sampleid"

#echo "$(date) 3. Start bam to fq: $sampleid"
#/NAS/wg_hl/Tools/samtools-1.9/samtools1.9/bin/samtools fastq -@ 20  $inputbam >  $clean_fq
#echo "$(date) 3. Finish bam to fq : $sampleid"

echo "$(date) 3. Start mapping: $sampleid"
$minimap2  --MD -ax map-ont -L -t 40 $ref $clean_fq > $sam
echo "$(date) 3. Finish mapping: $sampleid"

echo "$(date) 3. Start bam: $sampleid"
$samtools view -@ 20 -bS $sam -o $bam
echo "$(date) 3. Finish bam: $sampleid"

echo "$(date) 3. Start Delete sam: $sampleid"
rm $sam
echo "$(date) 3. Finish Delete sam: $sampleid"

echo "$(date) 3. Start sort: $sampleid"
$samtools sort  -@ 20 $bam -o $bam_sort
echo "$(date) 3. Finish sort: $sampleid"

echo "$(date) 4. Start add tag: $sampleid"
($samtools view -H "$bam_sort"; $samtools view "$bam_sort"|awk 'BEGIN{OFS="\t"}{$1="blood"$1; print $0}')|$samtools view -S -b - -o "$bam_sort_tag"
echo "$(date) 4. Finish add tag: $sampleid"

echo "$(date) 5. Start Delete bam: $sampleid"
rm $bam
rm $bam_sort
echo "$(date) 5. Finish Delete bam: $sampleid"

