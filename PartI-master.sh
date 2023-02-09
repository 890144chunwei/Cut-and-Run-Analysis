#!/usr/bin/bash
#SBATCH --mem 32768
#SBATCH -p short
#SBATCH --output=/storage/goodell/projects/chunweic/slurm_out/220726_peaks_ckit_wt_%j.out
#SBATCH -e /storage/goodell/projects/chunweic/slurm_out/220726_peaks_ckit_wt_%j.err # Standard output and error log
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=chunweic@bcm.edu # Email to which notifications will be sent

pwd; hostname; date

FASTQDIR="/storage/goodell/projects/chunweic/220601_Chunwei_CaR_molm13/220601"
OUTDIR="/storage/goodell/projects/chunweic/220630_Chunwei_CaR_LSK/Align"
HOMEDIR="/storage/goodell/home/chunweic"

trimmomatic PE -threads 1 -phred33 \
$FASTQDIR/WT_ckit_H3K27me3_S2_L001_R1_001.fastq $FASTQDIR/WT_ckit_H3K27me3_S2_L001_R2_001.fastq \
$OUTDIR/WT_ckit_H3K27me3_R1_paired.fastq $OUTDIR/WT_ckit_H3K27me3_R1_unpaired.fastq $OUTDIR/WT_ckit_H3K27me3_R2_paired.fastq $OUTDIR/WT_ckit_H3K27me3_R2_unpaired.fastq \
ILLUMINACLIP:$HOMEDIR/adapter_trim/TruSeq3-PE-2.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25
bowtie2 --dovetail --phred33 -x $HOMEDIR/mm10/BowtieIndex/mm10/mm10 -1 $OUTDIR/WT_ckit_H3K27me3_R1_paired.fastq -2 $OUTDIR/WT_ckit_H3K27me3_R2_paired.fastq -S $OUTDIR/WT_ckit_H3K27me3.sam
samtools view -S -b $OUTDIR/WT_ckit_H3K27me3.sam > $OUTDIR/WT_ckit_H3K27me3.bam
samtools sort $OUTDIR/WT_ckit_H3K27me3.bam $OUTDIR/WT_ckit_H3K27me3_sort
samtools index $OUTDIR/WT_ckit_H3K27me3_sort.bam

samtools rmdup -S $OUTDIR/WT_ckit_mIgG_sort.bam $OUTDIR/WT_ckit_mIgG_drm.bam
samtools sort $OUTDIR/WT_ckit_mIgG_drm.bam $OUTDIR/WT_ckit_mIgG_drm_sort
samtools index $OUTDIR/WT_ckit_mIgG_drm_sort.bam
bedtools bamtobed -i $OUTDIR/WT_ckit_mIgG_drm_sort.bam > $OUTDIR/WT_ckit_mIgG_H2AX1.bed
bedtools intersect -v -a $OUTDIR/WT_ckit_mIgG.bed -b $HOMEDIR/mm10-blacklist.v2.bed > $OUTDIR/WT_ckit_mIgG_bl.bed
bamCoverage -b $OUTDIR/WT_ckit_mIgG_drm_sort.bam -bl $HOMEDIR/mm10-blacklist.v2.bed -o $OUTDIR/WT_ckit_mIgG_cov.bw --normalizeUsing RPKM

Mark=("H2AZ" "H2AX" "SRCAP")
for target in ${Mark[@]}; do
  macs3 callpeak -g mm -f BAMPE -t $OUTDIR/Mut_ckit_H2AX1_drm_sort.bam $OUTDIR/Mut_ckit_H2AX2_drm_sort.bam -c $OUTDIR/Mut_ckit_mIgG_drm_sort.bam -q 0.05 -n $OUTDIR/Mut_ckit_H2AX --nomodel --keep-dup=all --call-summits
  macs3 callpeak -g mm -f BAMPE -t $OUTDIR/Mut_ckit_H2AZ1_drm_sort.bam $OUTDIR/Mut_ckit_H2AZ2_drm_sort.bam -c $OUTDIR/Mut_ckit_mIgG_drm_sort.bam -q 0.05 -n $OUTDIR/Mut_ckit_H2AZ --nomodel --keep-dup=all --call-summits
  macs3 callpeak -g mm -f BAMPE -t $OUTDIR/Mut_ckit_SRCAP1_drm_sort.bam $OUTDIR/Mut_ckit_SRCAP2_drm_sort.bam -c $OUTDIR/Mut_ckit_mIgG_drm_sort.bam -q 0.05 -n $OUTDIR/Mut_ckit_SRCAP --nomodel --keep-dup=all --call-summits
  annotatePeaks.pl $OUTDIR/Mut_ckit_H2AX_summits.bed mm10 -gtf $HOMEDIR/mm10/gencode.vM10.chr_patch_hapl_scaff.annotation.gtf > $OUTDIR/Mut_ckit_H2AX.txt
  annotatePeaks.pl $OUTDIR/Mut_ckit_H2AZ_summits.bed mm10 -gtf $HOMEDIR/mm10/gencode.vM10.chr_patch_hapl_scaff.annotation.gtf > $OUTDIR/Mut_ckit_H2AZ.txt
  annotatePeaks.pl $OUTDIR/Mut_ckit_SRCAP_summits.bed mm10 -gtf $HOMEDIR/mm10/gencode.vM10.chr_patch_hapl_scaff.annotation.gtf > $OUTDIR/Mut_ckit_SRCAP.txt
  cut -f 1-3 $OUTDIR/Mut_ckit_H2AX_peaks.narrowPeak > $OUTDIR/Mut_ckit_H2AX_peaks.bed
  cut -f 1-3 $OUTDIR/Mut_ckit_H2AZ_peaks.narrowPeak > $OUTDIR/Mut_ckit_H2AZ_peaks.bed
  cut -f 1-3 $OUTDIR/Mut_ckit_SRCAP_peaks.narrowPeak > $OUTDIR/Mut_ckit_SRCAP_peaks.bed
done
cat $OUTDIR/WT_ckit_H2AX_peaks.bed $OUTDIR/Mut_ckit_H2AX_peaks.bed > $OUTDIR/H2AX_ckit_combined_peaks.bed
bedtools sort -i $OUTDIR/H2AX_ckit_combined_peaks.bed > $OUTDIR/H2AX_ckit_combined_peaks_sort.bed
bedtools merge -i $OUTDIR/H2AX_ckit_combined_peaks_sort.bed > $OUTDIR/H2AX_ckit_combined_peaks_merge.bed

cat $OUTDIR/WT_ckit_H2AZ_peaks.bed $OUTDIR/Mut_ckit_H2AZ_peaks.bed > $OUTDIR/H2AZ_ckit_combined_peaks.bed
bedtools sort -i $OUTDIR/H2AZ_ckit_combined_peaks.bed > $OUTDIR/H2AZ_ckit_combined_peaks_sort.bed
bedtools merge -i $OUTDIR/H2AZ_ckit_combined_peaks_sort.bed > $OUTDIR/H2AZ_ckit_combined_peaks_merge.bed

cat $OUTDIR/WT_ckit_SRCAP_peaks.bed $OUTDIR/Mut_ckit_SRCAP_peaks.bed > $OUTDIR/SRCAP_ckit_combined_peaks.bed
bedtools sort -i $OUTDIR/SRCAP_ckit_combined_peaks.bed > $OUTDIR/SRCAP_ckit_combined_peaks_sort.bed
bedtools merge -i $OUTDIR/SRCAP_ckit_combined_peaks_sort.bed > $OUTDIR/SRCAP_ckit_combined_peaks_merge.bed

awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' $OUTDIR/H2AX_ckit_combined_peaks_merge.bed > $OUTDIR/H2AX_ckit_combined_merge.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' $OUTDIR/H2AZ_ckit_combined_peaks_merge.bed > $OUTDIR/H2AZ_ckit_combined_merge.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' $OUTDIR/SRCAP_ckit_combined_peaks_merge.bed > $OUTDIR/SRCAP_ckit_combined_merge.saf

featureCounts -T 7 -p -F SAF -a $OUTDIR/H2AX_ckit_combined_merge.saf -o $OUTDIR/WT_ckit_H2AX1_fc.txt $OUTDIR/WT_ckit_H2AX1_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/H2AX_ckit_combined_merge.saf -o $OUTDIR/WT_ckit_H2AX2_fc.txt $OUTDIR/WT_ckit_H2AX2_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/H2AX_ckit_combined_merge.saf -o $OUTDIR/Mut_ckit_H2AX1_fc.txt $OUTDIR/Mut_ckit_H2AX1_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/H2AX_ckit_combined_merge.saf -o $OUTDIR/Mut_ckit_H2AX2_fc.txt $OUTDIR/Mut_ckit_H2AX2_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/H2AZ_ckit_combined_merge.saf -o $OUTDIR/WT_ckit_H2AZ1_fc.txt $OUTDIR/WT_ckit_H2AZ1_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/H2AZ_ckit_combined_merge.saf -o $OUTDIR/WT_ckit_H2AZ2_fc.txt $OUTDIR/WT_ckit_H2AZ2_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/H2AZ_ckit_combined_merge.saf -o $OUTDIR/Mut_ckit_H2AZ1_fc.txt $OUTDIR/Mut_ckit_H2AZ1_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/H2AZ_ckit_combined_merge.saf -o $OUTDIR/Mut_ckit_H2AZ2_fc.txt $OUTDIR/Mut_ckit_H2AZ2_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/SRCAP_ckit_combined_merge.saf -o $OUTDIR/WT_ckit_SRCAP1_fc.txt $OUTDIR/WT_ckit_SRCAP1_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/SRCAP_ckit_combined_merge.saf -o $OUTDIR/WT_ckit_SRCAP2_fc.txt $OUTDIR/WT_ckit_SRCAP2_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/SRCAP_ckit_combined_merge.saf -o $OUTDIR/Mut_ckit_SRCAP1_fc.txt $OUTDIR/Mut_ckit_SRCAP1_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/SRCAP_ckit_combined_merge.saf -o $OUTDIR/Mut_ckit_SRCAP2_fc.txt $OUTDIR/Mut_ckit_SRCAP2_drm_sort.bam


