# Nextflow command used to run the pipeline

Using nextflow 21.10.3
Tissues
```bash
nextflow \
    run \
    nf-core/rnaseq \
        --input tissues_samplesheet.csv \
        --fasta Phygri_1.0.fasta \
        --gff Phygri1.3_gene_models_becker.gff \
        --save_reference \
        -profile biohpc_gen \
        --skip_stringtie \
        -process.maxRetries 2 \
        --outdir tissues_results \
        -r 3.9
```

Elicitors
```bash
nextflow \
    run \
    nf-core/rnaseq \
        --input elicitors_samplesheet.csv \
        --fasta Phygri_1.0.fasta \
        --gff Phygri1.3_gene_models_becker.gff \
        -profile biohpc_gen \
        --skip_stringtie \
        -process.maxRetries 2 \
        --outdir elicitors_results \
        -r 3.9
```