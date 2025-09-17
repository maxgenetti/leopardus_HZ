import os
import numpy as np

os.makedirs("./bash_jobs", exist_ok=True)
os.makedirs("./results", exist_ok=True)

demos = ['wedge','wedgeI','modern','wedgeS','wedgeIS','modernS']
sims = 10000
garbage = 5
ne0 = 0.10 #guttulus population density in individuals per square kilometer

seeds = np.random.choice(np.arange(1000, 2**31, dtype=np.int64), size=sims*len(demos), replace=False)

j = 1 #models are recorded as 1-indexed array
print("Models\tne0\tne1\tmig")
for i in list(range(0,sims)) :
    
    r=np.random.beta(a=2, b=3)
    ratio = 0.1 + r * (5.0 - 0.1)
    ne1 = round(ratio*ne0,4)

    r = np.random.beta(a=2, b=2)
    m = round(0.01 + r * (0.1111 - 0.01),4)  #1/9 would be complete replacement, 0.01 is 1 individual for guttulus at ratio of 0.8

    print(f"{j}-{j+5}\t{ne0}\t{ne1}\t{m}")
    for d in demos:
        seed = seeds[j-1] #corrects for j being 1-indexed for models
        out = open(f"./bash_jobs/{j}.sh", "w")
        out.write(f"python leoSELAM.py --seed {seed} --garbage {garbage} --ne0 {ne0} --ne1 {ne1} --mig {m} --demo {d} --model {j} 1>./results/model.{j}.tsv 2>./results/model.{j}.err\n")
        out.close()
        j+=1
