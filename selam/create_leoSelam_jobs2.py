import numpy as np
import random

demos = ['wedge','wedgeI','modern','wedgeS','wedgeIS','modernS']
sims = 10000
garbage = 5
ne0 = 0.10 #guttulus population density in individuals per square kilometer

ne1_r = np.arange(0.1,5.0000001,0.0001) #population density of geoffroyi indviduals relative to guttulus
m_r =   np.arange(0.01,0.111101,0.0001) #1/9 would be complete replacement, 0.01 is 1 individual for guttulus at ratio of 0.8
seeds = np.random.choice(np.arange(1000, 2**31, dtype=np.int64), size=sims*len(demos), replace=False)

j = 1 #models are recorded as 1-indexed array
print("Models\tne1\tne2\tmig")
for i in list(range(0,sims)) :
    ne1 = str(round(np.random.choice(ne1_r)*ne0,4))
    m = str(round(np.random.choice(m_r),4))
    print(f"{j}-{j+2}\t{ne0}\t{ne1}\t{m}")
    seed = seeds[j-1] #corrects for j being 1-indexed for models
    for d in demos:
        out = open(f"./bash_jobs/{j}.sh", "w")
        out.write(f"python leoSELAM.py --seeds {seed} --garbage {garbage} --ne0 {ne0} --ne1 {ne1} --mig {m} --demo {d} --model {j} 1>./results/model.{j}.tsv 2>./results/model.{j}.err\n")
        out.close()
        j+=1
