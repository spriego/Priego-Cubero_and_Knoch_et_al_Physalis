#  Methylseq

```bash
module load nextflow/23.10.0
```

Prepare genome

```bash
# I downloaded the lambda fasta seq from here https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1?report=fasta

cat /Phygri_1.0.fasta /phage_lambda.fasta > Phygri_1.0_and_lambda_phage.fasta

# Because the scaffold >tig00003611_subseq_846389:851390 Phygri_1.0 Physalis_grisea G20 gave some problems running the pipeline, I will remove it

# linearize the FASTA with awk, pipe to grep to filter for items of interest named in patterns.txt, then pipe to tr to delinearize":

awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' Phygri_1.0_and_lambda_phage.fasta | grep -v ">tig00003611_subseq_846389:851390 Phygri_1.0 Physalis_grisea G20" | tr "\t" "\n" > Phygri_1.0_lambda_modif.fasta 
```

### Run nextflow: 

Important, I called the samples for instance `root1` instead of `root_1`. The pipeline does not handle well the sample name and otherwise will only process 2 samples called `root` and `shoot`. 

Also, I have run into this issue https://github.com/nf-core/methylseq/issues/352 and it does not seem to be fixed at the time I am running the pipeline

I will follow the solution posted by @oliviapetrillo:

```
// inside prepare_genome.nf
SAMTOOLS_FAIDX([:], ch_fasta)
```

instead of

```
SAMTOOLS_FAIDX([[:], ch_fasta])
```

And then modifying SAMTOOLS_FAIDX:main.nf:
```
input: 
val(meta)
path(fasta)
```
instead of

```
input: 
tuple val(meta), path(fasta)
```

After changing the code above in the pipeline: 
```bash
nextflow \
	run \
	nf-core/methylseq \
	--input samplesheet.csv \
	--fasta Phygri_1.0_lambda_modif.fasta \
	--aligner bwameth \
	--comprehensive \
	--clip_r1 3 \
	--clip_r2 3 \
	--three_prime_clip_r1 3 \
	--three_prime_clip_r2 3 \
	--min_depth 5 \
	--outdir results \
	-profile biohpc_gen \
	-r 2.6.0
```


# MethylScore 

### Prepare combined bedGraph per sample 

Input for MethylScore:

```bash
# root
cat results/methyldackel/root1.markdup.sorted_C*.bedGraph > root_1.sorted.markDups_allC.bedGraph
cat results/methyldackel/root2.markdup.sorted_C*.bedGraph > root_2.sorted.markDups_allC.bedGraph
cat results/methyldackel/root5.markdup.sorted_C*.bedGraph > root_5.sorted.markDups_allC.bedGraph
# shoot
cat results/methyldackel/shoot1.markdup.sorted_C*.bedGraph > shoot_1.sorted.markDups_allC.bedGraph
cat results/methyldackel/shoot2.markdup.sorted_C*.bedGraph > shoot_2.sorted.markDups_allC.bedGraph
cat results/methyldackel/shoot3.markdup.sorted_C*.bedGraph > shoot_3.sorted.markDups_allC.bedGraph

```

Samplesheet:

```
root_1	root_1.sorted.markDups_allC.bedGraph
root_2	root_2.sorted.markDups_allC.bedGraph
root_5	root_5.sorted.markDups_allC.bedGraph
shoot_1	shoot_1.sorted.markDups_allC.bedGraph
shoot_2	shoot_2.sorted.markDups_allC.bedGraph
shoot_3	shoot_3.sorted.markDups_allC.bedGraph
```

### Run MethylScore:

```bash
nextflow \
	run \
	Computomics/MethylScore \
	--BEDGRAPH \
	--SAMPLE_SHEET=/samplesheet \
	--GENOME=Phygri_1.0_lambda_modif.fasta \
	-profile biohpc_gen
```


# Non-conversion rate

Calculated for all the contexts

```bash
for meth_file in *_allC.bedGraph
do
    filename=$(basename "$meth_file")
    output_file="ncr/${filename}.ncr"

	awk '$1 ~ /NC_001416.1/ {{m += $5}; {u += $6}} END {t = (m+u); print (m/t)*100}' $meth_file | awk -v filename="$filename" '{print $0, filename}' > $output_file
done

cat /dss/lxclscratch/08/ra57rut/physalis_wgbs/03_post_processing/a_ncr/*.ncr > all_samples_ncr.tsv
rm /dss/lxclscratch/08/ra57rut/physalis_wgbs/03_post_processing/a_ncr/*.ncr 
```

# Sup. Fig. 9B: Relative frequency of methylation rates

Based on Schultz et al 2012  Site methylation level = m / t

 m = The number of alignments/pairs reporting methylated bases
 u = The number of alignments/pairs reporting unmethylated bases
 t = Total of pairs reporting methylated and unmethylated bases

consider only sites covered more than 10

```bash
# file from MethylSeq Results 

for file in results/methyldackel/*.bedGraph;
do
	filename=$(basename "$file")
    output_file="rfmr/${file}.rfmr"

	grep -vE 'NC_001416.1|chrM|chrC' $file | awk -v 'NR>1 { if ( ($5+$6) > 9) print $5/($5+$6) }' OFMT='%.2g' | sort | uniq -c | awk -v OFS='\t' -v filename="$filename" '{print $0, filename}' > $output_file
done

cat rfmr/*.rfmr > rfmr/all_samples_relative_frequency_methylation_rates.tsv
rm rfmr/*.rfmr
```


# Mean methylation per window 

MethylDackel outputs the methylation percentage rounded to an integer on the 4th column

### For figure 4B
```bash
# 1. Get chromosomes length 
## samtools 1.20


samtools view -H results/bwameth/deduplicated/root1.markdup.sorted.bam|grep @SQ|sed 's/@SQ\tSN:\|LN://g' > Physalis_chromosomes_lenght.txt

# 2. Make windows
## bedtools 2.30.0
bedtools \
	makewindows \
	-g Physalis_chromosomes_lenght.txt \
	-w 100  > windows/Physalis_100bp_windows.bed

# 3. Mean methylation per window 
for file in results/methyldackel/*.bedGraph; do
    
	filename=$(basename "$file")
    output_file="windows/${filename}.100bp_window.mm"
	
    bedtools \
	map \
	-a windows/Physalis_100bp_windows.bed \
	-b "$file" \
	-g "Physalis_chromosomes_lenght.txt" \
	-c 4 \
	-o mean > $output_file
done	
```

### Sup. Fig. 8

Methylation per window for whole-genome representation (300 Kb)

```bash
# 1. Make windows
## bedtools 2.30.0
bedtools \
	makewindows \
	-g Physalis_chromosomes_lenght.txt \
	-w 300000  > windows/Physalis_300kb_windows.bed

# 2. Mean methylation per window 
for file in results/methyldackel/*.bedGraph; do
    
	filename=$(basename "$file")
    output_file="windows/${filename}.300kb_window.mm.tsv"
	
    bedtools \
	map \
	-a Physalis_300kb_windows.bed \
	-b "$file" \
	-g "Physalis_chromosomes_lenght.txt" \
	-c 4 \
	-o mean > $output_file
done	
```


### Sup. Fig. 9A: PCA on DMRs

```bash 
# MethylScore Results

for file in results/05DMRs/all/*.bed; do

	filename=$(basename "$file")
    output_file="dmrs/${filename}.position.bed"
	
    awk -F"\t" '{print $1 "\t" $2 "\t" $3}' $file > $output_file
done

# CG 
for file in results/methyldackel/*CpG.bedGraph; do

    filename=$(basename "$file")
    output_file="dmrs/${filename}.DMRs_Mean_Methylation.tsv"

	bedtools \
	map \
	-a dmrs/DMRs.CG.position.bed \
	-b "$file" \
	-g Physalis_chromosomes_lenght.txt \
	-c 4 \
	-o mean > $output_file
done

# CHG 
for file in results/methyldackel/*CHG.bedGraph; do

    filename=$(basename "$file")
    output_file="dmrs/${filename}.DMRs_Mean_Methylation.tsv"

	bedtools \
	map \
	-a dmrs/DMRs.CHG.position.bed \
	-b "$file" \
	-g Physalis_chromosomes_lenght.txt \
	-c 4 \
	-o mean > $output_file
done

# CHH 
for file in results/methyldackel/*CHH.bedGraph; do

    filename=$(basename "$file")
    output_file="dmrs/${filename}.DMRs_Mean_Methylation.tsv"

	bedtools \
	map \
	-a dmrs/DMRs.CHH.position.bed \
	-b "$file" \
	-g Physalis_chromosomes_lenght.txt \
	-c 4 \
	-o mean > $output_file
done
```


### Sup. Fig. 9B: Cs coverage in the cluster region 

```bash
# Make desired window

bedtools \
	makewindows \
	-g Physalis_chromosomes_lenght.txt \
	-w 1000  > windows/Physalis_2kb_windows.bed

# Filter for BGC
grep -w 'chr1' windows/Physalis_2kb_windows.bed | awk -v OFS='\t' '{ if(($2 >= 114700000) && ($3 <= 115120000)) { print } }' > windows/Physalis_2kb_windows_BGC.bed
```

```bash
# Single C coverage in the BGC locus
for file in results/methyldackel/*.bedGraph; do
    
	filename=$(basename "$file")
    output_file="coverage/${filename}.Cs_coverage_BGC.bed"
	
	grep -w 'chr1' $file | awk -v OFS='\t' '{ if(($2 >= 114700000) && ($3 <= 115120000)) { print $1, $2, $3, $5+$6 } }' > $output_file
done

for file in coverage/*.bed; do
    
	filename=$(basename "$file")
    output_file="coverage/${filename}.2kb_window.bed"
	
    bedtools \
	map \
	-a windows/Physalis_2kb_windows_BGC.bed \
	-b "$file" \
	-g "Physalis_chromosomes_lenght.txt" \
	-c 4 \
	-o mean > $output_file
done	
```


