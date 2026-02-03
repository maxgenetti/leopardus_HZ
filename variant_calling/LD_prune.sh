#!/bin/bash
# ld_prune.sh
# LD prune independently within each parental species

VCF_DIR="path/to/filtered_vcfs"
OUT_DIR="path/to/pruned_vcfs"
WINDOW=100  # kb
LD_THRESH=0.9
THREADS=8
CHROMS="A1 B1 C1 D1 E1 A2 B2 C2 D2 E2 A3 B3 C3 D3 E3 B4 D4"

# Split by parental species
for chrom in $CHROMS; do
    bcftools view $VCF_DIR/chr${chrom}.filt.1.vcf.gz \
        -S guttulus.samples.txt \
        -Oz -o $OUT_DIR/${chrom}.gut.vcf.gz
    bcftools index $OUT_DIR/${chrom}.gut.vcf.gz
    
    bcftools view $VCF_DIR/chr${chrom}.filt.1.vcf.gz \
        -S geoffroyi.samples.txt \
        -Oz -o $OUT_DIR/${chrom}.geo.vcf.gz
    bcftools index $OUT_DIR/${chrom}.geo.vcf.gz
done

# LD prune each species independently
for species in gut geo; do
    for chrom in $CHROMS; do
        bcftools +prune $OUT_DIR/${chrom}.${species}.vcf.gz \
            -e 'F_MISSING>=0.5' \
            -l $LD_THRESH \
            -w ${WINDOW}kb \
            -Oz -o $OUT_DIR/${chrom}.${species}.pruned.vcf.gz
        bcftools index $OUT_DIR/${chrom}.${species}.pruned.vcf.gz
    done
    
    # Concatenate chromosomes
    bcftools concat \
        $OUT_DIR/{A1,B1,C1,D1,E1,A2,B2,C2,D2,E2,A3,B3,C3,D3,E3,B4,D4}.${species}.pruned.vcf.gz \
        -Oz -o $OUT_DIR/${species}.pruned.vcf.gz
    bcftools index $OUT_DIR/${species}.pruned.vcf.gz
done

# Find shared sites between species (intersection)
bcftools isec \
    $OUT_DIR/geo.pruned.vcf.gz \
    $OUT_DIR/gut.pruned.vcf.gz \
    -p $OUT_DIR/isec_output

# Extract shared sites
bcftools query -f '%CHROM\t%POS\n' \
    $OUT_DIR/isec_output/0003.vcf \
    > $OUT_DIR/shared_pruned_sites.txt

# Subset original VCFs to shared pruned sites
for chrom in $CHROMS; do
    bcftools view $VCF_DIR/chr${chrom}.filt.1.vcf.gz \
        -R $OUT_DIR/shared_pruned_sites.txt \
        -Oz -o $OUT_DIR/${chrom}.final.vcf.gz
    bcftools index $OUT_DIR/${chrom}.final.vcf.gz
done

# Final concatenated VCF
bcftools concat \
    $OUT_DIR/{A1,B1,C1,D1,E1,A2,B2,C2,D2,E2,A3,B3,C3,D3,E3,B4,D4}.final.vcf.gz \
    -Oz -o $OUT_DIR/cohort.ld_pruned.vcf.gz
bcftools index $OUT_DIR/cohort.ld_pruned.vcf.gz