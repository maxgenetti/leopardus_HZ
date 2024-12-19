import numpy as np
import random
ne0_r = np.arange(0.005,0.1000001,0.0001)
ne1_r = np.arange(0.005,0.2000001,0.0001)
m_r =   np.arange(0.005,0.1000001,0.0001)

j = 1
print("Models\tne1\tne2\tmig")
for i in list(range(0,5000)) :
    ne0 = str(round(np.random.choice(ne0_r),4))
    ne1 = str(round(np.random.choice(ne1_r),4))
    m = str(round(np.random.choice(m_r),4))
    print(f"{j}-{j+5}\t{ne0}\t{ne1}\t{m}")
    for d in ['wedgeS','wedgeIS','modernS']:
        out = open(f"./bash_jobs/{j}.sh", "w")
        out.write(f"python leoSELAM.py --ne0 {ne0} --ne1 {ne1} --mig {m} --demo {d} --model {j} 1>./results/model.{j}.tsv 2>./results/model.{j}.err\n")
        out.close()
        j+=1
