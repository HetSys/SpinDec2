import matplotlib.pyplot as plt
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





f = open("output.txt", "r")
vals = f.read()

vals = vals.split("\n")
valsf = []
count = 0
all = []
for i in Ns:
    Nl = []
    for j in K:
        Kl = []
        for s in stabs:
            Sl = []
            for m in Ms:
                Ml = []
                for b in Bs:
                    Ml.append(float(vals[count]))
                    count = count +1
                Sl.append(Ml)
            Kl.append(Sl)
        Nl.append(Kl)
    all.append(Nl)

all = np.array(all)


plt.plot(all[:,0,0,0,0])
plt.show()
