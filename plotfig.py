import matplotlib.pyplot as plt 
import matplotlib.gridspec as grid 
import numpy as np 
import obspy 
data2 = []
file_ = open("build/out4.txt", "r") 
for idx, line in enumerate(file_.readlines()):
    sline = line.split(",")
    d = []
    for i in range(7):
        d.append(float(sline[i*2])+float(sline[i*2+1])*1j);
    data2.append(d)
    #if idx > 30 and idx < 60:
    #    print(idx, d[:3])
file_.close()

data2 = np.array(data2) 

data = data2

fig = plt.figure(1, figsize=(18, 18), dpi=300)

gs = grid.GridSpec(6, 1)
idxs = [4, 1, 5, 3, 0, 2]
for i in range(6):
    ax = fig.add_subplot((gs[i, 0]))
    d = np.fft.ifft(data2[:, idxs[i]])
    t = np.arange(len(d)) * 0.01
    ax.plot(t, np.real(d), alpha=0.5)
#plt.xlim((5, 15))
plt.savefig("img2.pdf")