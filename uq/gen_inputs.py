import numpy as np

f = open("uq_in.txt","r")
vals = f.read()
f.close()
valf = []
for i in vals.split("\n")[:-2]:
    valf.append(float(i.split("=")[-1]))
prob = vals.split("\n")[-2].split("=")[-1]


Ns = np.linspace(int(valf[0]),int(valf[1]),int(valf[2]),dtype=int)
K = np.linspace(valf[3],valf[4],int(valf[5]))
stabs = np.linspace(valf[6],valf[7],int(valf[8]))
Ms = np.linspace(valf[9],valf[10],int(valf[11]))
Bs = np.linspace(valf[12],valf[13],int(valf[14]))

count = 0

for i in Ns:
    for j in K:
        for s in stabs:
            for m in Ms:
                for b in Bs:
                    f = open("input_"+str(count)+".txt", "w")
                    f.write("""Concentration_min = 0.1
Concentration_max = 0.9
Domain_x_size = """+str(i)+"""
Domain_y_size = """+str(i)+"""
Mobility_A = """+str(m)+"""
Mobility_B = """+str(m)+"""
free_energy_gradient_parameter = """+str(j)+"""
Bulk_free_energy = """+str(b)+"""
Exitation_A = 0.1
Exitation_B = 0.2
Temperature_min = 800
Temperature_max = 1000
Problem = """+prob+"""
Stabilization_Term = """+str(s)+"""
Random_seed = -1
Checkpointing_interval = 5000
Max_time = 1e-1
time_step = 1e-4
dF_tolerance = 2.0
Use_input = 0
""")
                    f.close()
                    count = count + 1
