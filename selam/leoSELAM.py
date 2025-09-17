from collections import defaultdict
import numpy as np
import time
import argparse
import pandas as pd
from scipy.stats import ks_2samp
import sys
import os
import subprocess
import pickle
import random

def parseARGS () :
    parser = argparse.ArgumentParser()
    parser.add_argument('--demo', type=str, default="modern", action='store',help='demographic map')
    parser.add_argument('--output', type=str, default="out.tsv", action='store',help='output format file name in ./inputs/, default=out.tsv')
    parser.add_argument('--model', type=int, default=-1, action='store',help='model number')
    parser.add_argument('--ne0', type=float, default = 0.1, action='store', help='pop0 density km^-1')
    parser.add_argument('--ne1',type=float, default = 0.1, action='store', help='pop1 density km^-1')
    parser.add_argument('--mig', type=float, default = 0.1, action='store', help='dispersal into neighboring population')
    parser.add_argument('--seed', type=int, default=np.random.randint(1000, 2**31, dtype=np.int64), action='store',help='seed for SELAM')
    parser.add_argument('--garbage', type=int, default=5, action='store',help='garbage for SELAM, recommend 1-20 to manage memory usage')
    return parser.parse_args()

def read_chroms(chromFile = "./inputs/chroms.tsv"):
    chromDict = {}
    with open(chromFile) as f :
        for line in f :
            chrom = line.strip().split()[0]
            if chrom != 'X' :
                chromDict[chrom] = int(line.strip().split()[1])
    chromSizes = np.array(list(chromDict.values()))/50000000 #conver to morgans 2cM/Mb
    return np.round(chromSizes, 3)


def create_demo(demo, ne0, ne1, m, model) :
    df = pd.read_csv(f"./inputs/{demo}.txt",sep = "\t")
    area = df.loc[df['definition'] == 'size', ['1']].min().iloc[0]
    df['ne'] = df['species'].map({"geo":ne1,"gut":ne0,"hyb":ne0})
    df.loc[(df['definition'] == 'migration'), ['1']] = (df['1'] * m)
    df.loc[(df['definition'] == 'origin'), ['1']] = (df['1'] * m).round(5) ### added for edges
    df.loc[(df['definition'] == 'size') & (df['species'] == 'geo'), ['0']] = (df['0'] * ne1 * df['suitability']).round(0)
    df.loc[(df['definition'] == 'size') & (df['species'] == 'geo'), ['1']] = (df['1'] * ne1 * df['suitability']).round(0)
    df.loc[(df['definition'] == 'size') & (df['species'] == 'gut'), ['0']] = (df['0'] * ne0 * df['suitability']).round(0)
    df.loc[(df['definition'] == 'size') & (df['species'] == 'gut'), ['1']] = (df['1'] * ne0 * df['suitability']).round(0)
    df.loc[(df['definition'] == 'size') & (df['species'] == 'hyb'), ['0']] = (df['0'] * ne0 * df['suitability']).round(0)
    df.loc[(df['definition'] == 'size') & (df['species'] == 'hyb'), ['1']] = (df['1'] * ne0 * df['suitability']).round(0)
    filter_df = df[df['definition'] == 'migration'].copy()
    filter_df['area'] = area
    filter_df['filter'] = filter_df['suitability'] * filter_df['1'] * filter_df['ne'] * area < 1

    df = df.drop(df.loc[filter_df[filter_df['filter']==True].index].index)
    min_popSize = {}
    demo_out = open(f"./format_files/demography{model}.txt", "w")
    demo_out.write('pop0\tpop1\tsex\t0\t1\n')
    for index, row in df.iterrows():
        if row['definition'] == 'size':
            s1 = max(5,int(row['0'])) #max(4, int(row['0'])) ###min popsize of 4
            s2 = s1 ###min popsize of 4
            demo_out.write(f"{row['pop0']}\t{row['pop1']}\t{row['sex']}\t{s1}\t{s2}\n")
            min_popSize[row['pop0']] = min(5,int(s1/2)-1) ### set output as 1/4 of popsize to prevent it from breaking
        elif row['definition'] == 'origin':
            demo_out.write(f"{row['pop0']}\t{row['pop1']}\t{row['sex']}\t{row['0']:.0f}\t{row['1']:.3f}\n")
        elif row['definition'] == 'migration':
            recipient = df[(df['pop0'].astype(str)==str(row['pop0']))&(df['pop1'].astype(str)==str(row['pop0']))]['0'].values[0]
            donor =     df[(df['pop0'].astype(str)==str(row['pop1']))&(df['pop1'].astype(str)==str(row['pop1']))]['0'].values[0]
            if recipient > donor :
                mig = row['1'] * (donor/recipient)
            else :
                mig = row['1']

            demo_out.write(f"{row['pop0']}\t{row['pop1']}\t{row['sex']}\t0\t{mig:.5f}\n")

    demo_out.close()
    return(min_popSize)

def create_fifo(fifo_path):
     if not os.path.exists(fifo_path):
         os.mkfifo(fifo_path)
     else :
         print(f"FAILED TO MAKE {fifo_path}", file = sys.stderr)
         quit()

def create_output(fifo_path, model, min_popSize, output):
    popFile = pd.read_csv(f"./inputs/{output}", sep ="\t", header = None)
    popFile[2] = popFile[1].map(min_popSize)
    popFile[3] = popFile[1].map(min_popSize)
    popFile[4] = fifo_path
    popFile.to_csv(f"./format_files/output{model}.txt", sep="\t",index=False, header=None)
    return(popFile)

def create_viterbi(popFile):
    with open('./inputs/data.pkl', 'rb') as file:
        viterbi = pickle.load(file) #popData[pop][year]["1,1"] = list for each individual

    final_gen = popFile.iloc[popFile.index[-1]][0]
    final_pop = popFile.iloc[popFile.index[-1]][1]
    if popFile.iloc[popFile.index[-1]][3] != 0 :
        final_sex = 1
        final_ind = popFile.iloc[popFile.index[-1]][3]
    else :
        final_sex = 0
        final_ind = popFile.iloc[popFile.index[-1]][2]
    exit_state = [str(final_gen),str(final_pop),str(final_sex),str(final_ind)]
    return viterbi, exit_state

def k_s_test(ancestry_data, gen, pop, viterbi,summary_df) :
    for tract, sim in ancestry_data[gen][pop].items() :
        aa = sum(ancestry_data[gen][pop]['2,0'])
        ab = sum(ancestry_data[gen][pop]['1,1'])
        bb = sum(ancestry_data[gen][pop]['0,2'])
        ancestry_proportion = (aa+ab/2)/(aa+ab+bb)
        try :
            for year, anc_dict in viterbi[int(pop)].items() :
                for indv, lst in enumerate(anc_dict[tract]) :

                    ks_stat, p_value = ks_2samp(sim, lst)

                    #lst2 = random.choice(random.choice(list(random.choice(list(viterbi.values())).values()))[tract]) #randomly fit data to sim
                    #ks_stat2, p_value2 = ks_2samp(sim, lst2)
                    ks_stat2 = len(lst)
                    p_value2 = len(sim)

                    summary_df.loc[len(summary_df)] = [int(gen),pop,int(year),indv,tract,ancestry_proportion,ks_stat,p_value,ks_stat2,p_value2]

        except :
            #print(f"FAILED TO TEST\t-\tGEN:{gen}\tPOP:{pop}\tANC:{tract}\t-\tSIM_Len:{len(sim)}", file = sys.stderr)
            #return(pop)
            summary_df.loc[len(summary_df)] = [int(gen),pop,0,0,'0,2',ancestry_proportion,-1,-1,-1,-1]
            summary_df.loc[len(summary_df)] = [int(gen),pop,0,0,'1,1',ancestry_proportion,-1,-1,-1,-1]
            summary_df.loc[len(summary_df)] = [int(gen),pop,0,0,'2,0',ancestry_proportion,-1,-1,-1,-1]

    #if there is no admixture at 50 generations, then stop because the results should not change
    if int(gen)==55 and len(summary_df[(summary_df['GEN']==50)&(summary_df['SUM']>0)&(summary_df['SUM']<1)]) == 0 :
        return(gen)

def process_model_file(filepath, viterbi, exit_state, selam_process, model) :
    # Initialize data structures
    current_pop = -1
    current_gen = -1
    sex = -1
    ind = -1
    counts = np.array([0,0,0]) ### check for missing ancestries
    ancestry_data = {}
    maternal = []
    results_df = pd.DataFrame(columns = ["GEN","POP","YEAR","INDV","ANC","SUM","KST","P_v","KST2","P_v2"])
    with open(filepath, 'r') as fifo:
        while True:
            line = fifo.readline()
            parts = line.strip().split()
            if len(parts) < 6 :
                continue #empty line?
            elif parts[5] == "mtdna" : # end of individual. fills in missing anc tracks as 0. if last, process and exit
                maternal = []
                for idx, count in enumerate(counts) :
                    if count == 0 :
                        ancestry_data[generation][pop][["2,0","1,1","0,2"][idx]] += [0]
                counts = np.array([0,0,0])
                if exit_state == [generation, pop, sex, individual] :
                    x = k_s_test(ancestry_data, generation, pop, viterbi,results_df) #last pop and gen in outfifo
                    summarize_results(results_df, model)
                    print("DONE", file = sys.stderr)
                    return
                else :
                    continue #reset after paternal chromosome
            elif len(parts) != 9 : #
                if parts[6] == '1' : #gap after paternal chromosome, reset maternal
                    maternal = []
                continue #skip transition between parentals.

            generation, pop, sex, individual, chromosome, parent, ancestry, start, stop = parts
            if parent == "0" :
                maternal.append([int(ancestry), float(start), float(stop)])
            else :
                if generation not in ancestry_data.keys() :
                    if current_pop != -1 and current_gen != -1 :
                        x =  k_s_test(ancestry_data, current_gen, current_pop, viterbi,results_df)
                        if x :
                            summarize_results(results_df, model)
                            print(f"Terminating SELAM process due to simulation failure past threshold: {x}", file=sys.stderr)
                            selam_process.terminate()  # Gracefully terminate the process
                            selam_process.wait()  # Ensure the process has terminated
                            print("DONE", file = sys.stderr)
                            return
                        sys.stdout.flush() #flush after each generation
                    ancestry_data[generation] = {pop:{"2,0":[],"1,1":[],"0,2":[]}}
                    current_pop = pop
                    current_gen = generation
                elif pop not in ancestry_data[generation].keys() :
                    #new pop, same generation, so analyze previous pop and clear data
                    x =  k_s_test(ancestry_data, current_gen, current_pop, viterbi,results_df)
                    if x :
                        summarize_results(results_df, model)
                        print(f"Terminating SELAM process due to simulation failure past threshold: {x}", file=sys.stderr)
                        selam_process.terminate()  # Gracefully terminate the process
                        selam_process.wait()  # Ensure the process has terminated
                        print("DONE", file = sys.stderr)
                        return
                    ancestry_data[generation][pop] = {"2,0":[],"1,1":[],"0,2":[]}
                    current_pop = pop
                    current_gen = generation
                paternal_tract = [int(ancestry), float(start), float(stop)]
                a, b, c = calculate_overlap(maternal,paternal_tract)
                counts = counts+np.array([len(a),len(b),len(c)])
                ancestry_data[generation][pop]["2,0"] += a
                ancestry_data[generation][pop]["1,1"] += b
                ancestry_data[generation][pop]["0,2"] += c


def summarize_results(results_df, model) :

    results_df['GEN2'] = results_df['GEN']+results_df['YEAR']
    results_df.to_parquet(f'./results/model.{model}.parquet')

    summary_df = results_df.groupby('GEN2').agg(
        mean_p_value=('P_v', 'mean'),
        median_p_value=('P_v', 'median'),
        std_p_value=('P_v', 'std'),
        mean_p_value2=('P_v2', 'mean'),
        median_p_value2=('P_v2', 'median'),
        std_p_value2=('P_v2', 'std'),
        mean_ks_distance=('KST', 'mean'),
        median_ks_distance=('KST', 'median'),
        std_KST=('KST', 'std'),
        max_ks_distance=('KST', 'max'),
        mean_ks_distance2=('KST2', 'mean'),
        median_ks_distance2=('KST2', 'median'),
        std_KST2=('KST2', 'std'),
        max_ks_distance2=('KST2', 'max'),
        count_subpopulations=('GEN2', 'count')  # Count how many subpopulations are used
    )
    summary_df.to_csv(sys.stdout, sep='\t', index_label='GEN2')


def calculate_overlap(mother_tracts, father_tract):
    father_state, f_start, f_stop = father_tract
    overlap_lengths = {
        0 : [],
        1 : [],
        2 : []
    }
    for mother_state, m_start, m_stop in mother_tracts:
        # Calculate the overlap range
        start = max(f_start, m_start)
        stop = min(f_stop, m_stop)

        if start < stop:  # Check if there is an actual overlap
            # Update the corresponding state combination with the overlap length
            overlap_lengths[mother_state+father_state].append(stop - start)

    return overlap_lengths[0],overlap_lengths[1],overlap_lengths[2]

def main():
    start_time = time.time()

    selam_path = "~/bin/SELAM/src/SELAM"

    args = parseARGS()
    ne0 = args.ne0
    ne1 = args.ne1
    m = args.mig
    model = args.model
    demo = args.demo
    output = args.output
    seed = args.seed
    garbage = args.garbage
    print(f"{model}\t{ne0}\t{ne1}\t{m}",file = sys.stderr) #store model parameters

    ####set up directory
    if not os.path.exists(f"./results/") :
        os.makedirs(f"./results/")
    if not os.path.exists(f"./output_fifos/") :
        os.makedirs(f"./output_fifos/")
    if not os.path.exists(f"./format_files/") :
        os.makedirs(f"./format_files/")

    #####make output fifo
    min_popSize = create_demo(demo, ne0, ne1, m, model)
    chromSizes = read_chroms()
    fifo_path = f"./output_fifos/output_fifo_{model}.txt"
    create_fifo(fifo_path)
    popFile = create_output(fifo_path, model, min_popSize, output)
    viterbi, exit_state = create_viterbi(popFile)

    ##########run command
    pChrom = ' '.join([str(x) for x in chromSizes])
    command = f"{selam_path} --seed {seed} --garbage {garbage} -d ./format_files/demography{model}.txt -o ./format_files/output{model}.txt -c {len(chromSizes)+1} {pChrom} 0 >>selam_runtime.txt"
    print(f"Starting : {command}", file = sys.stderr)

    selam_process = subprocess.Popen(command, shell=True, preexec_fn=os.setsid)
    process_model_file(fifo_path, viterbi, exit_state, selam_process, model)

    ##########clean up directories
    os.popen(f"rm {fifo_path}")
    os.popen(f"rm ./format_files/demography{model}.txt")
    os.popen(f"rm ./format_files/output{model}.txt")

    print(f"Finished in {time.time()-start_time} seconds", file = sys.stderr)

if __name__ == "__main__":
    main()
