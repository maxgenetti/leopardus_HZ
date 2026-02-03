#!/usr/bin/env python3
"""
Convert VCF to ancestry HMM input format
Reads LD-pruned, biallelic SNP VCF and converts to aHMM panel format
"""

import gzip
import sys
import argparse
import numpy as np
from collections import defaultdict

class VCF2aHMM:
    """
    Convert filtered VCF to ancestry HMM format
    
    Args:
        vcf: Path to VCF file (filtered for LD, indels, biallelic SNPs)
        pop: Population assignment file (parent1, parent2, admixed)
        rate: Mean recombination rate in cM/Mb
        dist: Minimum distance between sites in bp
        minDif: Minimum allele frequency difference between parents
        minGT: Minimum genotype calls per parent population
        minDP: Minimum proportion of admixed samples with coverage
        geno: Use genotypes instead of read counts for admixed samples
    
    Returns:
        Ancestry HMM input panel (written to stdout)
    """

    def __init__(self):
        args = self.get_args()
        self.vName = args.vcf
        self.pName = args.pop
        self.rate = args.rate * 0.00000001  # Convert cM/Mb to morgan/bp
        self.minFreq = args.minDif
        self.minGT = args.minGT
        self.minDP = args.minDP
        self.dist = args.dist
        self.geno = args.geno
        
        # VCF parsing attributes
        self.format_idx = int()
        self.ad_idx = int()
        self.gt_idx = int()
        
        # Output tracking
        self.siteFreq = []
        self.current_pos = 0
        self.panel_lines = []
        
        # Sample assignments
        self.ancestry = {}
        self.header = []
        self.chrom = None

    def get_args(self):
        parser = argparse.ArgumentParser(
            description='Convert VCF to ancestry HMM input format'
        )
        parser.add_argument('--vcf', type=str, required=True,
                           help='Path to VCF file')
        parser.add_argument('--pop', type=str, required=True,
                           help='Population assignment file')
        parser.add_argument('--rate', type=float, default=2.0,
                           help='Mean recombination rate (cM/Mb)')
        parser.add_argument('--dist', type=int, default=1000,
                           help='Minimum distance between sites (bp)')
        parser.add_argument('--minDif', type=float, default=0.2,
                           help='Minimum allele frequency difference')
        parser.add_argument('--minGT', type=int, default=20,
                           help='Minimum genotype calls per parent')
        parser.add_argument('--minDP', type=float, default=0.5,
                           help='Minimum admixed sample coverage')
        parser.add_argument('-geno', action='store_true', default=False,
                           help='Use genotypes instead of read counts')
        return parser.parse_args()

    def read_samples(self):
        """Parse population assignment file"""
        self.ancestry['parent1'] = []
        self.ancestry['parent2'] = []
        self.ancestry['admixed'] = []
        
        with open(self.pName) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                    
                pop = parts[0]
                samples = parts[1].split(',')
                
                for sample in samples:
                    try:
                        # Handle numeric indices
                        idx = int(sample) + self.format_idx
                        self.ancestry[pop].append(idx)
                    except ValueError:
                        # Handle sample names
                        if sample in self.header:
                            self.ancestry[pop].append(self.header.index(sample))

    def read_vcf(self):
        """Read and process VCF file"""
        # Open VCF
        if self.vName.endswith("gz"):
            vcf = gzip.open(self.vName, "rt")
        else:
            vcf = open(self.vName, "r")
        
        # Parse header
        for line in vcf:
            if line.startswith('#CHROM'):
                self.header = line.strip().split()
                self.format_idx = self.header.index('FORMAT')
                break
        
        # Read sample assignments
        self.read_samples()
        
        # Get chromosome from first data line
        first_line = vcf.readline()
        fields = first_line.strip().split()
        self.chrom = fields[0]
        
        # Get FORMAT field indices
        format_fields = fields[self.format_idx].split(':')
        self.gt_idx = format_fields.index('GT')
        
        if not self.geno:
            if 'AD' not in format_fields:
                print("ERROR: No AD field in FORMAT, use -geno option", file=sys.stderr)
                sys.exit(1)
            self.ad_idx = format_fields.index('AD')
        
        # Process first line
        self.process_site(fields)
        
        # Process remaining lines
        for line in vcf:
            fields = line.strip().split()
            self.process_site(fields)
        
        vcf.close()

    def process_site(self, fields):
        """Process a single VCF site"""
        pos = int(fields[1])
        
        # Check minimum distance
        if pos - self.current_pos <= self.dist:
            return
        
        # Calculate genetic distance
        morgans = (pos - self.current_pos) * self.rate
        
        # Count alleles in each parent population
        parent1_counts = [0, 0]
        parent2_counts = [0, 0]
        
        for idx in self.ancestry['parent1']:
            gt = fields[idx].split(':')[self.gt_idx]
            try:
                parent1_counts[int(gt[0])] += 1
                parent1_counts[int(gt[-1])] += 1
            except (ValueError, IndexError):
                continue
        
        if sum(parent1_counts) < self.minGT:
            return
        
        for idx in self.ancestry['parent2']:
            gt = fields[idx].split(':')[self.gt_idx]
            try:
                parent2_counts[int(gt[0])] += 1
                parent2_counts[int(gt[-1])] += 1
            except (ValueError, IndexError):
                continue
        
        if sum(parent2_counts) < self.minGT:
            return
        
        # Calculate allele frequencies
        p1_freq = parent1_counts[0] / sum(parent1_counts)
        p2_freq = parent2_counts[0] / sum(parent2_counts)
        freq_diff = abs(p1_freq - p2_freq)
        
        # Check minimum frequency difference
        if freq_diff < self.minFreq:
            return
        
        # Get admixed sample data
        admixed_data = []
        
        if self.geno:
            # Use genotypes
            for idx in self.ancestry['admixed']:
                gt = fields[idx].split(':')[self.gt_idx]
                alleles = [0, 0]
                try:
                    alleles[int(gt[0])] += 1
                    alleles[int(gt[-1])] += 1
                except (ValueError, IndexError):
                    pass
                admixed_data.extend(alleles)
        else:
            # Use allele depths
            for idx in self.ancestry['admixed']:
                ad = fields[idx].split(':')[self.ad_idx].split(',')
                alleles = [0, 0]
                try:
                    alleles[0] = int(ad[0])
                    alleles[1] = int(ad[1])
                except (ValueError, IndexError):
                    pass
                admixed_data.extend(alleles)
        
        # Check coverage in admixed samples
        sample_coverage = [admixed_data[i] + admixed_data[i+1] 
                          for i in range(0, len(admixed_data), 2)]
        coverage_rate = np.count_nonzero(sample_coverage) / len(self.ancestry['admixed'])
        
        if coverage_rate < self.minDP:
            return
        
        # Write output line
        output = [
            self.chrom,
            str(pos),
            str(parent1_counts[0]),
            str(parent1_counts[1]),
            str(parent2_counts[0]),
            str(parent2_counts[1]),
            f"{morgans:.4E}"
        ] + [str(x) for x in admixed_data]
        
        self.panel_lines.append('\t'.join(output))
        self.siteFreq.append(freq_diff)
        self.current_pos = pos

    def write_output(self):
        """Write panel to stdout and stats to stderr"""
        # Write panel
        for line in self.panel_lines:
            print(line)
        
        # Write statistics
        print(f"Sites in Panel             : {len(self.siteFreq)}", file=sys.stderr)
        if self.siteFreq:
            print(f"Mean Frequency Difference  : {np.mean(self.siteFreq):.4f}", file=sys.stderr)
            print(f"Frequency Difference stDev : {np.std(self.siteFreq):.4f}", file=sys.stderr)


def main():
    converter = VCF2aHMM()
    converter.read_vcf()
    converter.write_output()


if __name__ == "__main__":
    main()