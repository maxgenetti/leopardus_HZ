#!/bin/bash
# call_variants.sh
# Joint variant calling with GATK HaplotypeCaller

REFERENCE="path/to/reference.fa"
BAM_DIR="path/to/bams"
OUT_DIR="path/to/vcfs"
THREADS=8

# Define chromosomes
declare -A chroms
chroms[59326]="A1"
chroms[59327]="B1"
chroms[59328]="C1"
chroms[59329]="D1"
chroms[59330]="E1"
chroms[59331]="A2"
chroms[59332]="B2"
chroms[59333]="C2"
chroms[59334]="D2"
chroms[59335]="E2"
chroms[59336]="A3"
chroms[59337]="B3"
chroms[59338]="C3"
chroms[59339]="D3"
chroms[59340]="E3"
chroms[59341]="B4"
chroms[59342]="D4"
chroms[59343]="X"

# Define chromosome lengths
declare -A chrom_lengths
chrom_lengths[59326]=239694388
chrom_lengths[59327]=205836458
chrom_lengths[59328]=222052948
chrom_lengths[59329]=115783437
chrom_lengths[59330]=61591894
chrom_lengths[59331]=169709481
chrom_lengths[59332]=152606360
chrom_lengths[59333]=158996348
chrom_lengths[59334]=88472993
chrom_lengths[59335]=62031975
chrom_lengths[59336]=140630183
chrom_lengths[59337]=148130213
chrom_lengths[59338]=153706148
chrom_lengths[59339]=94884351
chrom_lengths[59340]=42047855
chrom_lengths[59341]=142838203
chrom_lengths[59342]=94461142
chrom_lengths[59343]=128930408
chrom_lengths[28320]=16738

BLOCK_SIZE=10000000  # 10MB

# Process autosomes (skip MT for now)
for chrom_id in {59326..59343}; do
    chrom=${chroms[$chrom_id]}
    length=${chrom_lengths[$chrom_id]}
    
    echo "Processing chromosome ${chrom} (NC_0${chrom_id}.1), length: ${length}"
    
    # Calculate number of blocks
    num_blocks=$(( (length + BLOCK_SIZE - 1) / BLOCK_SIZE ))
    
    # Call variants for each 10MB block
    for block in $(seq 0 $((num_blocks - 1))); do
        start=$((block * BLOCK_SIZE + 1))
        end=$(( (block + 1) * BLOCK_SIZE ))
        if [ $end -gt $length ]; then
            end=$length
        fi
        
        echo "  Block ${block}: ${start}-${end}"
        
        gatk HaplotypeCaller \
            -R $REFERENCE \
            -I $BAM_DIR/*.bam \
            -L NC_0${chrom_id}.1:${start}-${end} \
            -O $OUT_DIR/chr${chrom}.block${block}.vcf.gz \
            --native-pair-hmm-threads $THREADS
    done
    
    # Combine all blocks for this chromosome
    bcftools concat \
        $OUT_DIR/chr${chrom}.block*.vcf.gz \
        -Oz -o $OUT_DIR/chr${chrom}.raw.vcf.gz
    
    bcftools index $OUT_DIR/chr${chrom}.raw.vcf.gz
    
    # Clean up block files
    rm $OUT_DIR/chr${chrom}.block*.vcf.gz*
done