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

lins = [Ns,K,stabs,Ms,Bs]
lins_size = [len(Ns),len(K),len(stabs),len(Ms),len(Bs)]
dims = ["N","kappa","Stabilisation term","Mobility","Bulk free energy"]
dimp = dims[np.argmax(lins_size)]

to_plot = lins[np.argmax(lins_size)]

f = open("output.txt", "r")
vals = f.read()

vals = vals.split("\n")
valsf = []
count = 0
all = []
all2 = []
for i in Ns:
    Nl = []
    Nl2 = []
    for j in K:
        Kl = []
        Kl2 = []
        for s in stabs:
            Sl = []
            Sl2 = []
            for m in Ms:
                Ml = []
                Ml2 = []
                for b in Bs:
                    if prob[2:-1] == "Both":
                        Ml.append(float(vals[count]))
                        count = count +1
                        Ml2.append(float(vals[count]))
                        count = count +1
                    else:
                        Ml.append(float(vals[count]))
                        count = count +1
                Sl.append(Ml)
                Sl2.append(Ml2)
            Kl.append(Sl)
            Kl2.append(Sl2)
        Nl.append(Kl)
        Nl2.append(Kl2)
    all.append(Nl)
    all2.append(Nl2)

all = np.array(all)
all2 = np.array(all2)




if prob[2:-1] == "Both":
    plt.plot(to_plot,np.squeeze(all),label = "Constant")
    plt.plot(to_plot,np.squeeze(all2),label = "Spectral")
else:
    plt.plot(to_plot,np.squeeze(all),label = prob[2:-1])
plt.legend()
plt.xlabel(dimp)
plt.ylabel("Critical time step (s)")
plt.show()
