# 1. Run nextflow pipeline

```bash
nextflow \
	run \
	nf-core/atacseq \
	--input samplesheet.csv \
	--fasta Phygri_1.0.fasta.gz \
	--gff Phygri1.3_gene_models_becker.gff \
	--aligner bowtie2 \
	--read_length 150 \
	--with_control \
	--blacklist physalis_blacklist.bed \
	--mito_name 'chrM_1|chrM_2|chrM_3|chrM_4|chrM_5|chrM_6|chrM_7|chrC' \
	--ataqv_mito_reference 'chrM_7' \
	-profile biohpc_gen \
	-process.maxRetries 2 \
	--outdir results \
	-r 2.1.2
```    

# 2. Deeptools

Compute matrix
deeptools 3.5.5


```bash
computeMatrix scale-regions \
    -R physalis_becker_annotation/expression_quantiles/physalis_root_*.bed \
    -S results/bowtie2/merged_library/bigwig/ROOTS_*.bigWig results/bowtie2/merged_library/bigwig/INPUT_REP1.mLb.clN.bigWig \
    -b 3000 \
    -a 3000 \
    -p 12 \
    --missingDataAsZero \
    --regionBodyLength 3000 \
    -o deeptools/computematrix/root_scaled.gz \
    --outFileNameMatrix deeptools/computematrix/root_BW_scaled.tab \
    --outFileSortedRegions deeptools/computematrix/root_BW_scaled_genes.bed
```
The same code was used for leaf samples but using the expression quantiles of the 2-w-old shoot sample

```bash
plotHeatmap \
      -m deeptools/computematrix/root_scaled.gz \
      -out deeptools/plotheatmat/root_metaplots.pdf \
      --colorMap Oranges \
      --labelRotation 45 \
      --legendLocation upper-right \
      --averageTypeSummaryPlot median \
      --perGroup \
      --plotFileFormat pdf
```


# 3. Peak calling

Root samples

```bash
macs3 \
	hmmratac \
	-f BAMPE \
	-n root_all \
	-i results/bowtie2/merged_library/ROOTS_REP1.mLb.clN.sorted.bam results/bowtie2/merged_library/ROOTS_REP2.mLb.clN.sorted.bam results/bowtie2/merged_library/ROOTS_REP3.mLb.clN.sorted.bam results/bowtie2/merged_library/ROOTS_REP4.mLb.clN.sorted.bam \
	--outdir macs3/hmmratac/root
```

Leaf samples

```bash
macs3 \
	hmmratac \
	-f BAMPE \
	-n leaf_all \
	-i results/bowtie2/merged_library/LEAVES_REP1.mLb.clN.sorted.bam results/bowtie2/merged_library/LEAVES_REP2.mLb.clN.sorted.bam results/bowtie2/merged_library/LEAVES_REP3.mLb.clN.sorted.bam results/bowtie2/merged_library/LEAVES_REP4.mLb.clN.sorted.bam \
	--outdir macs3/hmmratac/leaf
```

Filter score > 100

```bash
awk '$5 > 100' macs3/hmmratac/leaf/leaf_all_accessible_regions.narrowPeak > macs3/hmmratac/leaf/leaf_all_accessible_regions.filtered.narrowPeak

awk '$5 > 100' macs3/hmmratac/root/root_all_accessible_regions.narrowPeak > macs3/hmmratac/root/root_all_accessible_regions.filtered.narrowPeak
```


### 3.1 Consensus peak file


Consensus file used on R to plot the PCA and the boxplot in Fig 2
```bash
bedtools intersect -a macs3/hmmratac/leaf/leaf_all_accessible_regions.filtered.narrowPeak -b macs3/hmmratac/root/root_all_accessible_regions.filtered.narrowPeak > macs3/bedtools/intersect_f.narrowPeak

bedtools intersect -a macs3/hmmratac/leaf/leaf_all_accessible_regions.filtered.narrowPeak -b macs3/hmmratac/root/root_all_accessible_regions.filtered.narrowPeak -v > macs3/bedtools/intersect_leaf_private.narrowPeak

bedtools intersect -a macs3/hmmratac/root/root_all_accessible_regions.filtered.narrowPeak -b macs3/hmmratac/leaf/leaf_all_accessible_regions.filtered.narrowPeak -v > macs3/bedtools/intersect_root_private.narrowPeak

cat macs3/bedtools/intersect_f.narrowPeak macs3/bedtools/intersect_leaf_private.narrowPeak /macs3/bedtools/intersect_root_private.narrowPeak | sort -k1,1 -k2,2n > macs3/bedtools/leaf_and_root_concensus.narrowPeak 
```

### 3.2 Get raw counts on the concensus peaks per sample

```bash
for bam_file in results/bowtie2/merged_library/*.bam
do
    filename=$(basename "$bam_file")
    output_file="macs3/bedtools/${filename}.concensus.narrowPeak.counts"

    bedtools \
        coverage \
        -a macs3/bedtools/leaf_and_root_concensus.narrowPeak \
        -b $bam_file \
        -counts > $output_file
done
```

# 4. Homer 

## 4.1 Annotate peaks
Annotate peaks of leaf and root sample detected with `hmmratac`

Same code applied to `root_all_accessible_regions.filtered.narrowPeak`

```bash
annotatePeaks.pl \
    macs3/hmmratac/leaf/leaf_all_accessible_regions.filtered.narrowPeak \
    Phygri_1.0.fasta.gz \
    -gtf Phygri1.3_gene_models.becker_agat_fixed.gtf > homer/annotate/leaf/leaf_all_accessible_regions.filtered.narrowPeak.annotation
```

## 4.2 Enrichment analysis


`root_sample_RESC_peaks.homer` = peaks detected in the root-expressed sub-cluster (RESC) in the file `root_all_accessible_regions.filtered.narrowPeak`

```bash
findMotifsGenome.pl \
    root_sample_RESC_peaks.homer \
    Phygri_1.0.fasta \
    homer/enrichment/root/root_all_results \
    -size 500 \
    -mask \
    -p 6
```

Same code was used for leaf samples. `leaf_sample_SESC_peaks.homer` = peaks detected in the shoot-expressed sub-cluster (SESC) in the file `leaf_all_accessible_regions.filtered.narrowPeak`