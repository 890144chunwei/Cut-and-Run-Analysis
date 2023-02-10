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

for FILE in $PROJECTDIR/*R1.fastq ;
do
  trimmomatic PE -threads 1 -phred33 $PROJECTDIR/${FILE%_R1.fastq}_R1.fastq $PROJECTDIR/${FILE%_R1.fastq}_R2.fastq \
  $OUTDIR/${FILE%_R1.fastq}_R1_paired.fastq $OUTDIR/${FILE%_R1.fastq}_R1_unpaired.fastq $OUTDIR/${FILE%_R1.fastq}_R2_paired.fastq $OUTDIR/${FILE%_R1.fastq}_R2_unpaired.fastq \
  ILLUMINACLIP:$HOMEDIR/adapter_trim/NexteraPE-PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25
  bowtie2 --dovetail --phred33 -x $HOMEDIR/mm10/BowtieIndex/mm10/mm10 -1 $OUTDIR/${FILE%_R1.fastq}_R1_paired.fastq -2 $OUTDIR/${FILE%_R1.fastq}_R2_paired.fastq -S $OUTDIR/${FILE%_R1.fastq}.sam
  samtools view -S -b $OUTDIR/${FILE%_R1.fastq}.sam > $OUTDIR/${FILE%_R1.fastq}.bam
  samtools rmdup -S $OUTDIR/${FILE%_R1.fastq}.bam $OUTDIR/${FILE%_R1.fastq}_drm.bam
  samtools sort $OUTDIR/${FILE%_R1.fastq}_drm.bam $OUTDIR/${FILE%_R1.fastq}_sort
  samtools index $OUTDIR/${FILE%_R1.fastq}_sort.bam
  bedtools bamtobed -i $OUTDIR/${FILE%_R1.fastq}_shift_sort.bam > $OUTDIR/${FILE%_R1.fastq}.bed
  bedtools intersect -v -a $OUTDIR/${FILE%_R1.fastq}.bed -b $HOMEDIR/mm10-blacklist.v2.bed > $OUTDIR/${FILE%_R1.fastq}_bl.bed
  bamCoverage -b $OUTDIR/${FILE%_R1.fastq}_shift_sort.bam -bl $HOMEDIR/mm10-blacklist.v2.bed -o $OUTDIR/${FILE%_R1.fastq}_cov.bw --normalizeUsing RPKM
done
rm $PROJECTDIR/*fastq $PROJECTDIR/*zip $PROJECTDIR/*shift.bam $PROJECTDIR/*drm.bam $OUTDIR/*unpaired.fastq

Mark=("H2AZ" "H2AX" "SRCAP")
for target in ${Mark[@]}; do
  macs3 callpeak -g mm -f BAMPE -t $OUTDIR/Mut_ckit_$target_drm_sort.bam -c $OUTDIR/Mut_ckit_mIgG_drm_sort.bam -q 0.05 -n $OUTDIR/Mut_ckit_$target --nomodel --keep-dup=all --call-summits
  annotatePeaks.pl $OUTDIR/Mut_ckit_$target_summits.bed mm10 -gtf $HOMEDIR/mm10/gencode.vM10.chr_patch_hapl_scaff.annotation.gtf > $OUTDIR/Mut_ckit_$target.txt
  cut -f 1-3 $OUTDIR/Mut_ckit_$target_peaks.narrowPeak > $OUTDIR/Mut_ckit_$target_peaks.bed
  cat $OUTDIR/WT_ckit_$target_peaks.bed $OUTDIR/Mut_ckit_$target_peaks.bed > $OUTDIR/$target_ckit_combined_peaks.bed
  bedtools sort -i $OUTDIR/H$target_ckit_combined_peaks.bed > $OUTDIR/$target_ckit_combined_peaks_sort.bed
  bedtools merge -i $OUTDIR/$target_ckit_combined_peaks_sort.bed > $OUTDIR/$target_ckit_combined_peaks_merge.bed
  awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' $OUTDIR/$target_ckit_combined_peaks_merge.bed > $OUTDIR/$target_ckit_combined_merge.saf
done

featureCounts -T 7 -p -F SAF -a $OUTDIR/SRCAP_ckit_combined_merge.saf -o $OUTDIR/WT_ckit_SRCAP1_fc.txt $OUTDIR/WT_ckit_SRCAP1_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/SRCAP_ckit_combined_merge.saf -o $OUTDIR/WT_ckit_SRCAP2_fc.txt $OUTDIR/WT_ckit_SRCAP2_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/SRCAP_ckit_combined_merge.saf -o $OUTDIR/Mut_ckit_SRCAP1_fc.txt $OUTDIR/Mut_ckit_SRCAP1_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/SRCAP_ckit_combined_merge.saf -o $OUTDIR/Mut_ckit_SRCAP2_fc.txt $OUTDIR/Mut_ckit_SRCAP2_drm_sort.bam
