#!/usr/bin/env python3
"""
Ancestry HMM workflow for chromosome-scale VCFs
Runs iterative ancestry inference with filtering of introgressed sites
"""

import numpy as np
import pandas as pd
import os
import sys
import subprocess
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Run ancestry HMM workflow on chromosome-scale VCFs')
    parser.add_argument('--chroms', '-c', type=str, required=True,
                        help='File containing list of chromosome names (one per line)')
    parser.add_argument('--directory', '-d', type=str, required=True,
                        help='Output directory name')
    parser.add_argument('--step', '-a', type=int, default=1,
                        help='Step in process to start at (1-7)')
    parser.add_argument('--threads', '-t', type=int, default=40,
                        help='Number of parallel threads')
    parser.add_argument('--samples', '-s', type=str, required=True,
                        help='File with list of admixed samples to analyze')
    parser.add_argument('--dist', type=int, default=5000,
                        help='Minimum distance between sites in bp')
    parser.add_argument('--diff', type=float, default=0.1,
                        help='Minimum allele frequency difference between parental populations')
    parser.add_argument('--minGT', type=int, default=30,
                        help='Minimum number of genotype calls in each parent population')
    parser.add_argument('--vcf-dir', type=str, default='/scratch/mgenetti/meowMix/vcf/pruned',
                        help='Directory containing LD-pruned VCFs')
    parser.add_argument('--parent1', type=str, default='geoffroyi.30.txt',
                        help='File with parent1 sample names')
    parser.add_argument('--parent2', type=str, default='guttulus.42.txt',
                        help='File with parent2 sample names')
    return parser.parse_args()


def read_sample_file(filename):
    """Read sample names from file, one per line"""
    samples = []
    with open(filename, 'r') as f:
        for line in f:
            samples.append(line.strip())
    return samples


def setup_directories(outdir, topdir, chroms):
    """Create output directory structure"""
    if not os.path.exists(f"{outdir}/{topdir}"):
        os.makedirs(f"{outdir}/{topdir}")
        os.makedirs(f"{outdir}/{topdir}/samples")
    
    for step in [1, 2, 3, 4, 5, 6, "V", "V2"]:
        step_dir = f"{outdir}/{topdir}/step{step}"
        if not os.path.exists(step_dir):
            os.makedirs(step_dir)
            # Create chromosome subdirectories for steps that need them
            if step in [1, 3, 5]:
                for chrom in chroms:
                    os.makedirs(f"{step_dir}/chr{chrom}")


def create_pop_files(sample, parent1, parent2, outdir, topdir):
    """Create population assignment files for a sample"""
    # Main sample population file
    pop_file = f"{outdir}/{topdir}/samples/{sample}.pop"
    with open(pop_file, 'w') as f:
        f.write(f"parent1\t{','.join([x for x in parent1 if x != sample])}\n")
        f.write(f"parent2\t{','.join([x for x in parent2 if x != sample])}\n")
        f.write(f"admixed\t{sample}")
    
    # Sample file for ancestry HMM
    sample_file = f"{outdir}/{topdir}/samples/{sample}.sample"
    with open(sample_file, 'w') as f:
        f.write(f"{sample}\t2")
    
    # If sample is a parent, create additional files for jackknife analysis
    if sample in parent1 + parent2:
        pop_file = f"{outdir}/{topdir}/samples/{sample}_sample.pop"
        with open(pop_file, 'w') as f:
            f.write(f"parent1\t{','.join(parent1)}\n")
            f.write(f"parent2\t{','.join(parent2)}\n")
            f.write(f"admixed\t{sample}")
        
        sample_file = f"{outdir}/{topdir}/samples/{sample}_sample.sample"
        with open(sample_file, 'w') as f:
            f.write(f"{sample}_sample\t2")


def step1_naive_run(args, outdir, samples, parent1, parent2, chroms):
    """Step 1: Run naive ancestry HMM without filtering parental introgression"""
    print(f"STARTING STEP 1", file=sys.stderr)
    
    topdir = args.directory
    step = 1
    
    # Create command files
    v2a_cmds = open(f"{outdir}/{topdir}/v2a{step}.txt", "w")
    v2a_fix_cmds = open(f"{outdir}/{topdir}/v2a{step+1}.txt", "w")
    ahmm_cmds = open(f"{outdir}/{topdir}/ahmm{step}.txt", "w")
    
    # Initial ancestry proportions
    p = 0.5
    q = 0.5
    
    for sample in samples:
        # Create population assignment files
        create_pop_files(sample, parent1, parent2, outdir, topdir)
        
        # Build VCF to ancestry HMM conversion commands for each chromosome
        panel_files = []
        for chrom in chroms:
            vcf_path = f"{args.vcf_dir}/{chrom}.final.vcf.gz"
            panel_out = f"{outdir}/{topdir}/step{step}/chr{chrom}/{sample}.panel"
            panel_err = f"{outdir}/{topdir}/step{step}/chr{chrom}/{sample}.panel.err"
            
            v2a_cmds.write(
                f"python vcf2aHMM.py "
                f"--vcf {vcf_path} "
                f"--dist {args.dist} "
                f"--pop {outdir}/{topdir}/samples/{sample}.pop "
                f"--rate 2 "
                f"--minGT {args.minGT} "
                f"--minDP 0.50 "
                f"--minDif {args.diff} "
                f"1>{panel_out} 2>{panel_err}\n"
            )
            panel_files.append(panel_out)
        
        # Concatenate chromosome panels
        all_panels = ' '.join(panel_files)
        final_panel = f"{outdir}/{topdir}/step{step}/{sample}.panel"
        v2a_fix_cmds.write(f"cat {all_panels} > {final_panel} && rm {all_panels}\n")
        
        # Run ancestry HMM
        ahmm_cmds.write(
            f"cd {outdir}/{topdir}/step{step}/ && "
            f"ancestry_hmm "
            f"-i {final_panel} "
            f"-s {outdir}/{topdir}/samples/{sample}.sample "
            f"-e 1e-3 --ne 10000 -a 2 {q} {p} "
            f"-p 0 10000 {q} -p 1 -100 {p} "
            f"1>{sample}.out 2>{sample}.err\n"
        )
    
    v2a_cmds.close()
    v2a_fix_cmds.close()
    ahmm_cmds.close()
    
    # Execute commands in parallel
    run_parallel(f"{outdir}/{topdir}/v2a{step}.txt", args.threads, f"{outdir}/{topdir}/log.v2a{step}")
    run_parallel(f"{outdir}/{topdir}/v2a{step+1}.txt", args.threads, f"{outdir}/{topdir}/log.v2a{step+1}")
    run_parallel(f"{outdir}/{topdir}/ahmm{step}.txt", args.threads, f"{outdir}/{topdir}/log.ahmm{step}")
    
    print(f"FINISHED STEP 1", file=sys.stderr)


def run_parallel(cmd_file, threads, logfile):
    """Run commands in parallel using GNU parallel"""
    command = f"parallel -j {threads} --joblog {logfile} < {cmd_file}"
    subprocess.run(command, shell=True, check=True)


def extract_introgressed_sites(posterior_file, sample, threshold=0.1):
    """
    Extract introgressed sites from ancestry HMM posterior probabilities
    Returns list of (chrom, pos) tuples for sites showing introgression
    """
    sites = []
    
    try:
        data = pd.read_csv(posterior_file, delimiter="\t")
        
        # Calculate introgression probability
        # Adjust based on actual ancestry_hmm posterior format
        for idx, row in data.iterrows():
            try:
                chrom = row['chr']
                pos = row['pos']
                
                # Calculate probability of non-native ancestry
                # This assumes columns like '2,0' for parent2 ancestry
                introgression_prob = 0.0
                
                # Sum non-native ancestry probabilities
                # Adjust column names based on actual output
                if '2,0' in row:
                    introgression_prob += row['2,0']
                if '0,2' in row:
                    introgression_prob += row['0,2']
                
                if introgression_prob > threshold:
                    sites.append((sample, chrom, pos))
            except:
                continue
    except Exception as e:
        print(f"Error processing {posterior_file}: {e}", file=sys.stderr)
    
    return sites


def main():
    args = get_args()
    
    # Set output directory
    outdir = "/scratch/mgenetti/meowMix/ahmm"
    
    # Read input files
    chroms = read_sample_file(args.chroms)
    admixed_samples = read_sample_file(args.samples)
    parent1 = read_sample_file(args.parent1)
    parent2 = read_sample_file(args.parent2)
    
    # Combine all samples (admixed + parents) and remove duplicates
    all_samples = sorted(list(set(admixed_samples + parent1 + parent2)))
    
    # Setup directory structure
    setup_directories(outdir, args.directory, chroms)
    
    # Run requested step
    step = args.step
    
    if step == 1:
        step1_naive_run(args, outdir, all_samples, parent1, parent2, chroms)
        step = 2
    
    if step == 2:
        step2_update_ancestry(args, outdir, all_samples, parent1, parent2, chroms)
        step = 3
    
    if step == 3:
        step3_filter_introgression(args, outdir, all_samples, parent1, parent2, chroms)
        step = 4
    
    if step == 4:
        step4_update_ancestry(args, outdir, all_samples, parent1, parent2, chroms)
        step = 5
    
    if step == 5:
        step5_filter_introgression(args, outdir, all_samples, parent1, parent2, chroms)
        step = 6
    
    if step == 6:
        step6_final_ancestry(args, outdir, all_samples, parent1, parent2, chroms)
        step = 7
    
    if step == 7:
        step7_summarize(args, outdir, all_samples, parent1, parent2)


def step2_update_ancestry(args, outdir, samples, parent1, parent2, chroms):
    """Step 2: Update ancestry proportions based on step 1 results"""
    if args.step == 2:
        # If starting at step 2, run the previous AHMM commands
        run_parallel(f"{outdir}/{args.directory}/ahmm1.txt", args.threads, 
                    f"{outdir}/{args.directory}/log.ahmm1")
    
    print(f"STARTING STEP 2", file=sys.stderr)
    
    topdir = args.directory
    step = 2
    
    ahmm_cmds = open(f"{outdir}/{topdir}/ahmm{step}.txt", "w")
    
    for sample in samples:
        try:
            # Read posterior probabilities from step 1
            data = pd.read_csv(f"{outdir}/{topdir}/step1/{sample}.posterior", delimiter="\t")
            data['post'] = data['2,0'] + 0.5 * data['1,1']
            mean = data["post"].mean()
        except:
            print(f"ERROR: Could not read {outdir}/{topdir}/step1/{sample}.posterior", file=sys.stderr)
            continue
        
        # Update ancestry proportions
        q = round(max(0.001, min(mean, 0.999)), 3)
        p = round(1.00 - q, 3)
        
        ahmm_cmds.write(
            f"cd {outdir}/{topdir}/step{step}/ && "
            f"ancestry_hmm "
            f"-i {outdir}/{topdir}/step1/{sample}.panel "
            f"-s {outdir}/{topdir}/samples/{sample}.sample "
            f"-e 1e-3 --ne 10000 -a 2 {q} {p} "
            f"-p 0 10000 {q} -p 1 -100 {p} "
            f"1>{sample}.out 2>{sample}.err\n"
        )
    
    ahmm_cmds.close()
    run_parallel(f"{outdir}/{topdir}/ahmm{step}.txt", args.threads, 
                f"{outdir}/{topdir}/log.ahmm{step}")
    
    print(f"FINISHED STEP 2", file=sys.stderr)


def step3_filter_introgression(args, outdir, samples, parent1, parent2, chroms):
    """Step 3: Filter introgressed sites in parental samples"""
    if args.step == 3:
        run_parallel(f"{outdir}/{args.directory}/ahmm2.txt", args.threads,
                    f"{outdir}/{args.directory}/log.ahmm2")
    
    print(f"STARTING STEP 3", file=sys.stderr)
    
    topdir = args.directory
    step = 3
    
    # Extract introgressed sites from all parental samples
    all_introgressed = []
    for sample in parent1 + parent2:
        posterior_file = f"{outdir}/{topdir}/step2/{sample}.posterior"
        sites = extract_introgressed_sites(posterior_file, sample, threshold=0.1)
        all_introgressed.extend(sites)
    
    # Write filter file
    filter_file = f"{outdir}/{topdir}/step{step}/all.filter"
    with open(filter_file, 'w') as f:
        for sample, chrom, pos in all_introgressed:
            f.write(f"{sample}\t{chrom}\t{pos}\n")
    
    # Re-run vcf2aHMM with filtered sites
    v2a_cmds = open(f"{outdir}/{topdir}/v2a{step}.txt", "w")
    v2a_fix_cmds = open(f"{outdir}/{topdir}/v2a{step+1}.txt", "w")
    ahmm_cmds = open(f"{outdir}/{topdir}/ahmm{step}.txt", "w")
    
    for sample in samples:
        # Get updated ancestry proportions
        try:
            data = pd.read_csv(f"{outdir}/{topdir}/step2/{sample}.posterior", delimiter="\t")
            data['post'] = data['2,0'] + 0.5 * data['1,1']
            mean = data["post"].mean()
        except:
            mean = 0.5
        
        q = round(max(0.001, min(mean, 0.999)), 3)
        p = round(1.00 - q, 3)
        
        # Build VCF conversion commands with filtering
        panel_files = []
        for chrom in chroms:
            vcf_path = f"{args.vcf_dir}/{chrom}.final.vcf.gz"
            panel_out = f"{outdir}/{topdir}/step{step}/chr{chrom}/{sample}.panel"
            panel_err = f"{outdir}/{topdir}/step{step}/chr{chrom}/{sample}.panel.err"
            
            v2a_cmds.write(
                f"python vcf2aHMM.py "
                f"--vcf {vcf_path} "
                f"--dist {args.dist} "
                f"--pop {outdir}/{topdir}/samples/{sample}.pop "
                f"--filter {outdir}/{topdir}/step{step}/all.filter "
                f"--rate 2 "
                f"--minGT {args.minGT} "
                f"--minDP 0.50 "
                f"--minDif {args.diff} "
                f"-fp "
                f"1>{panel_out} 2>{panel_err}\n"
            )
            panel_files.append(panel_out)
        
        # Concatenate chromosome panels
        all_panels = ' '.join(panel_files)
        final_panel = f"{outdir}/{topdir}/step{step}/{sample}.panel"
        v2a_fix_cmds.write(f"cat {all_panels} > {final_panel} && rm {all_panels}\n")
        
        # Run ancestry HMM
        ahmm_cmds.write(
            f"cd {outdir}/{topdir}/step{step}/ && "
            f"ancestry_hmm "
            f"-i {final_panel} "
            f"-s {outdir}/{topdir}/samples/{sample}.sample "
            f"-e 1e-3 --ne 10000 -a 2 {q} {p} "
            f"-p 0 10000 {q} -p 1 -100 {p} "
            f"1>{sample}.out 2>{sample}.err\n"
        )
    
    v2a_cmds.close()
    v2a_fix_cmds.close()
    ahmm_cmds.close()
    
    run_parallel(f"{outdir}/{topdir}/v2a{step}.txt", args.threads, f"{outdir}/{topdir}/log.v2a{step}")
    run_parallel(f"{outdir}/{topdir}/v2a{step+1}.txt", args.threads, f"{outdir}/{topdir}/log.v2a{step+1}")
    run_parallel(f"{outdir}/{topdir}/ahmm{step}.txt", args.threads, f"{outdir}/{topdir}/log.ahmm{step}")
    
    print(f"FINISHED STEP 3", file=sys.stderr)


def step4_update_ancestry(args, outdir, samples, parent1, parent2, chroms):
    """Step 4: Update ancestry proportions after filtering"""
    if args.step == 4:
        run_parallel(f"{outdir}/{args.directory}/ahmm3.txt", args.threads,
                    f"{outdir}/{args.directory}/log.ahmm3")
    
    print(f"STARTING STEP 4", file=sys.stderr)
    
    topdir = args.directory
    step = 4
    
    ahmm_cmds = open(f"{outdir}/{topdir}/ahmm{step}.txt", "w")
    
    for sample in samples:
        try:
            data = pd.read_csv(f"{outdir}/{topdir}/step3/{sample}.posterior", delimiter="\t")
            data['post'] = data['2,0'] + 0.5 * data['1,1']
            mean = data["post"].mean()
        except:
            print(f"ERROR: Could not read {outdir}/{topdir}/step3/{sample}.posterior", file=sys.stderr)
            continue
        
        q = round(max(0.001, min(mean, 0.999)), 3)
        p = round(1.00 - q, 3)
        
        ahmm_cmds.write(
            f"cd {outdir}/{topdir}/step{step}/ && "
            f"ancestry_hmm "
            f"-i {outdir}/{topdir}/step3/{sample}.panel "
            f"-s {outdir}/{topdir}/samples/{sample}.sample "
            f"-e 1e-3 --ne 10000 -a 2 {q} {p} "
            f"-p 0 10000 {q} -p 1 -100 {p} "
            f"1>{sample}.out 2>{sample}.err\n"
        )
    
    ahmm_cmds.close()
    run_parallel(f"{outdir}/{topdir}/ahmm{step}.txt", args.threads,
                f"{outdir}/{topdir}/log.ahmm{step}")
    
    print(f"FINISHED STEP 4", file=sys.stderr)


def step5_filter_introgression(args, outdir, samples, parent1, parent2, chroms):
    """Step 5: Second round of filtering introgressed sites"""
    if args.step == 5:
        run_parallel(f"{outdir}/{args.directory}/ahmm4.txt", args.threads,
                    f"{outdir}/{args.directory}/log.ahmm4")
    
    print(f"STARTING STEP 5", file=sys.stderr)
    
    topdir = args.directory
    step = 5
    
    # Extract introgressed sites
    all_introgressed = []
    for sample in parent1 + parent2:
        posterior_file = f"{outdir}/{topdir}/step4/{sample}.posterior"
        sites = extract_introgressed_sites(posterior_file, sample, threshold=0.1)
        all_introgressed.extend(sites)
    
    # Write filter file
    filter_file = f"{outdir}/{topdir}/step{step}/all.filter"
    with open(filter_file, 'w') as f:
        for sample, chrom, pos in all_introgressed:
            f.write(f"{sample}\t{chrom}\t{pos}\n")
    
    # Re-run with updated filters
    v2a_cmds = open(f"{outdir}/{topdir}/v2a{step}.txt", "w")
    v2a_fix_cmds = open(f"{outdir}/{topdir}/v2a{step+1}.txt", "w")
    ahmm_cmds = open(f"{outdir}/{topdir}/ahmm{step}.txt", "w")
    
    for sample in samples:
        try:
            data = pd.read_csv(f"{outdir}/{topdir}/step4/{sample}.posterior", delimiter="\t")
            data['post'] = data['2,0'] + 0.5 * data['1,1']
            mean = data["post"].mean()
        except:
            mean = 0.5
        
        q = round(max(0.001, min(mean, 0.999)), 3)
        p = round(1.00 - q, 3)
        
        panel_files = []
        for chrom in chroms:
            vcf_path = f"{args.vcf_dir}/{chrom}.final.vcf.gz"
            panel_out = f"{outdir}/{topdir}/step{step}/chr{chrom}/{sample}.panel"
            panel_err = f"{outdir}/{topdir}/step{step}/chr{chrom}/{sample}.panel.err"
            
            v2a_cmds.write(
                f"python vcf2aHMM.py "
                f"--vcf {vcf_path} "
                f"--dist {args.dist} "
                f"--pop {outdir}/{topdir}/samples/{sample}.pop "
                f"--filter {outdir}/{topdir}/step{step}/all.filter "
                f"--rate 2 "
                f"--minGT {args.minGT} "
                f"--minDP 0.50 "
                f"--minDif {args.diff} "
                f"-fp "
                f"1>{panel_out} 2>{panel_err}\n"
            )
            panel_files.append(panel_out)
        
        all_panels = ' '.join(panel_files)
        final_panel = f"{outdir}/{topdir}/step{step}/{sample}.panel"
        v2a_fix_cmds.write(f"cat {all_panels} > {final_panel} && rm {all_panels}\n")
        
        ahmm_cmds.write(
            f"cd {outdir}/{topdir}/step{step}/ && "
            f"ancestry_hmm "
            f"-i {final_panel} "
            f"-s {outdir}/{topdir}/samples/{sample}.sample "
            f"-e 1e-3 --ne 10000 -a 2 {q} {p} "
            f"-p 0 10000 {q} -p 1 -100 {p} "
            f"1>{sample}.out 2>{sample}.err\n"
        )
    
    v2a_cmds.close()
    v2a_fix_cmds.close()
    ahmm_cmds.close()
    
    run_parallel(f"{outdir}/{topdir}/v2a{step}.txt", args.threads, f"{outdir}/{topdir}/log.v2a{step}")
    run_parallel(f"{outdir}/{topdir}/v2a{step+1}.txt", args.threads, f"{outdir}/{topdir}/log.v2a{step+1}")
    run_parallel(f"{outdir}/{topdir}/ahmm{step}.txt", args.threads, f"{outdir}/{topdir}/log.ahmm{step}")
    
    print(f"FINISHED STEP 5", file=sys.stderr)


def step6_final_ancestry(args, outdir, samples, parent1, parent2, chroms):
    """Step 6: Final ancestry estimation and Viterbi decoding"""
    if args.step == 6:
        run_parallel(f"{outdir}/{args.directory}/ahmm5.txt", args.threads,
                    f"{outdir}/{args.directory}/log.ahmm5")
    
    print(f"STARTING STEP 6", file=sys.stderr)
    
    topdir = args.directory
    step = 6
    
    ahmm_cmds = open(f"{outdir}/{topdir}/ahmm{step}.txt", "w")
    ahmmv_cmds = open(f"{outdir}/{topdir}/ahmmV2.txt", "w")
    
    sum_df = pd.DataFrame(columns=["Sample", "Anc", "Gen", "Sites", "Distance", "Difference", "minGT", "color"])
    
    for sample in samples:
        if sample in parent2:
            color = "red"
        elif sample in parent1:
            color = "blue"
        else:
            color = "green"
        
        try:
            # Count sites
            result = subprocess.run(
                f"wc -l {outdir}/{topdir}/step5/{sample}.panel",
                shell=True, capture_output=True, check=True
            )
            sites = int(result.stdout.split()[0])
            
            # Read generation estimate
            with open(f"{outdir}/{topdir}/step5/{sample}.out", 'r') as f:
                lines = f.readlines()
                gen = float(lines[3].strip().split()[1])
            
            # Calculate mean ancestry
            data = pd.read_csv(f"{outdir}/{topdir}/step5/{sample}.posterior", delimiter="\t")
            data['post'] = data['2,0'] + 0.5 * data['1,1']
            mean = data["post"].mean()
        except Exception as e:
            print(f"ERROR processing {sample}: {e}", file=sys.stderr)
            continue
        
        q = round(max(0.001, min(mean, 0.999)), 3)
        p = round(1.00 - q, 3)
        
        # Viterbi decoding
        ahmmv_cmds.write(
            f"cd {outdir}/{topdir}/stepV/ && "
            f"ancestry_hmm -v "
            f"-i {outdir}/{topdir}/step5/{sample}.panel "
            f"-s {outdir}/{topdir}/samples/{sample}.sample "
            f"-e 1e-3 --ne 10000 -a 2 {q} {p} "
            f"-p 0 10000 {q} -p 1 -100 {p} "
            f"1>{sample}.out 2>{sample}.err\n"
        )
        
        # Final ancestry estimation
        ahmm_cmds.write(
            f"cd {outdir}/{topdir}/step{step}/ && "
            f"ancestry_hmm "
            f"-i {outdir}/{topdir}/step5/{sample}.panel "
            f"-s {outdir}/{topdir}/samples/{sample}.sample "
            f"-e 1e-3 --ne 10000 -a 2 {q} {p} "
            f"-p 0 10000 {q} -p 1 -100 {p} "
            f"1>{sample}.out 2>{sample}.err\n"
        )
        
        sum_df.loc[len(sum_df)] = [sample, mean, gen, sites, args.dist, args.diff, args.minGT, color]
    
    # Process parental samples separately
    for sample in parent1 + parent2:
        if sample in parent2:
            color = "orange"
        elif sample in parent1:
            color = "purple"
        else:
            color = "black"
        
        try:
            result = subprocess.run(
                f"wc -l {outdir}/{topdir}/step5/{sample}_sample.panel",
                shell=True, capture_output=True, check=True
            )
            sites = int(result.stdout.split()[0])
            
            with open(f"{outdir}/{topdir}/step5/{sample}_sample.out", 'r') as f:
                lines = f.readlines()
                gen = float(lines[3].strip().split()[1])
            
            data = pd.read_csv(f"{outdir}/{topdir}/step5/{sample}_sample.posterior", delimiter="\t")
            data['post'] = data['2,0'] + 0.5 * data['1,1']
            mean = data["post"].mean()
        except Exception as e:
            print(f"ERROR processing {sample}_sample: {e}", file=sys.stderr)
            continue
        
        q = round(max(0.001, min(mean, 0.999)), 3)
        p = round(1.00 - q, 3)
        
        ahmmv_cmds.write(
            f"cd {outdir}/{topdir}/stepV2/ && "
            f"ancestry_hmm -v "
            f"-i {outdir}/{topdir}/step5/{sample}_sample.panel "
            f"-s {outdir}/{topdir}/samples/{sample}_sample.sample "
            f"-e 1e-3 --ne 10000 -a 2 {q} {p} "
            f"-p 0 10000 {q} -p 1 -100 {p} "
            f"1>{sample}_sample.out 2>{sample}_sample.err\n"
        )
        
        ahmm_cmds.write(
            f"cd {outdir}/{topdir}/step{step}/ && "
            f"ancestry_hmm "
            f"-i {outdir}/{topdir}/step5/{sample}_sample.panel "
            f"-s {outdir}/{topdir}/samples/{sample}_sample.sample "
            f"-e 1e-3 --ne 10000 -a 2 {q} {p} "
            f"-p 0 10000 {q} -p 1 -100 {p} "
            f"1>{sample}_sample.out 2>{sample}_sample.err\n"
        )
        
        sum_df.loc[len(sum_df)] = [f"{sample}_sample", mean, gen, sites, args.dist, args.diff, args.minGT, color]
    
    ahmm_cmds.close()
    ahmmv_cmds.close()
    
    sum_df.to_csv(f'{outdir}/{topdir}/summary.step5.csv', index=False, sep="\t")
    
    run_parallel(f"{outdir}/{topdir}/ahmm{step}.txt", args.threads, f"{outdir}/{topdir}/log.ahmm{step}")
    run_parallel(f"{outdir}/{topdir}/ahmmV2.txt", args.threads, f"{outdir}/{topdir}/log.ahmmV2")
    
    print(f"FINISHED STEP 6", file=sys.stderr)


def step7_summarize(args, outdir, samples, parent1, parent2):
    """Step 7: Final summary statistics"""
    print(f"STARTING STEP 7", file=sys.stderr)
    
    topdir = args.directory
    step = 7
    
    sum_df = pd.DataFrame(columns=["Sample", "Anc", "Gen", "Sites", "Distance", "Difference", "minGT", "color"])
    
    for sample in samples:
        if sample in parent2:
            color = "red"
        elif sample in parent1:
            color = "blue"
        else:
            color = "green"
        
        try:
            result = subprocess.run(
                f"wc -l {outdir}/{topdir}/step5/{sample}.panel",
                shell=True, capture_output=True, check=True
            )
            sites = int(result.stdout.split()[0])
            
            with open(f"{outdir}/{topdir}/step6/{sample}.out", 'r') as f:
                lines = f.readlines()
                gen = float(lines[3].strip().split()[1])
            
            data = pd.read_csv(f"{outdir}/{topdir}/step6/{sample}.posterior", delimiter="\t")
            data['post'] = data['2,0'] + 0.5 * data['1,1']
            mean = data["post"].mean()
            
            sum_df.loc[len(sum_df)] = [sample, mean, gen, sites, args.dist, args.diff, args.minGT, color]
        except Exception as e:
            print(f"ERROR processing {sample}: {e}", file=sys.stderr)
    
    for sample in parent1 + parent2:
        if sample in parent2:
            color = "orange"
        elif sample in parent1:
            color = "purple"
        else:
            color = "black"
        
        try:
            result = subprocess.run(
                f"wc -l {outdir}/{topdir}/step5/{sample}_sample.panel",
                shell=True, capture_output=True, check=True
            )
            sites = int(result.stdout.split()[0])
            
            with open(f"{outdir}/{topdir}/step6/{sample}_sample.out", 'r') as f:
                lines = f.readlines()
                gen = float(lines[3].strip().split()[1])
            
            data = pd.read_csv(f"{outdir}/{topdir}/step6/{sample}_sample.posterior", delimiter="\t")
            data['post'] = data['2,0'] + 0.5 * data['1,1']
            mean = data["post"].mean()
            
            sum_df.loc[len(sum_df)] = [f"{sample}_sample", mean, gen, sites, args.dist, args.diff, args.minGT, color]
        except Exception as e:
            print(f"ERROR processing {sample}_sample: {e}", file=sys.stderr)
    
    sum_df.to_csv(f'{outdir}/{topdir}/summary.step6.csv', index=False, sep="\t")
    
    print(f"FINISHED STEP 7", file=sys.stderr)


if __name__ == "__main__":
    main()