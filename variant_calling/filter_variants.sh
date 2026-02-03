#!/bin/bash
# filter_variants.sh
# Filter VCFs for biallelic SNPs, MAF, depth

IN_DIR="path/to/raw_vcfs"
OUT_DIR="path/to/filtered_vcfs"
THREADS=8

for vcf in $IN_DIR/chr*.raw.vcf.gz; do
    chrom=$(basename $vcf .raw.vcf.gz | sed 's/chr//')
    
    bcftools filter --threads $THREADS -g 5 -Ou $vcf \
    | bcftools view --threads $THREADS \
        -i 'INFO/DP<=2000 && INFO/AN>=40' \
        -m2 -M2 \
        -q 0.01 -Q 0.99 \
        -v snps -Ou \
    | bcftools annotate --threads $THREADS \
        --rename-chrs rename.chroms.txt \
        -Oz -o $OUT_DIR/chr${chrom}.filt.vcf.gz
    
    bcftools index $OUT_DIR/chr${chrom}.filt.vcf.gz
    
    # Set missing genotypes (DP=0) to ./.
    bcftools +setGT $OUT_DIR/chr${chrom}.filt.vcf.gz \
        -t q -n . -i 'FMT/DP=0' \
        -Oz -o $OUT_DIR/chr${chrom}.filt.1.vcf.gz
    
        # Clean up intermediate files
    rm $OUT_DIR/chr${chrom}.filt.vcf.gz*

done