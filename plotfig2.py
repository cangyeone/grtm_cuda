import matplotlib.pyplot as plt 
import matplotlib.gridspec as grid 
import numpy as np 
import obspy 
data2 = []
file_ = open("build/signal.txt", "r") 
for idx, line in enumerate(file_.readlines()):
    sline = line.split(",")
    sline = [float(s) for s in sline[:-1]]
    data2.append(sline)
file_.close()

data2 = np.array(data2) 

data = data2

fig = plt.figure(1, figsize=(36, 18), dpi=300)

gs = grid.GridSpec(6, 2)
idxs = [4, 1, 5, 3, 0, 2]
for i in range(6):
    ax = fig.add_subplot((gs[i, 0]))

    #tr = obspy.read(f"/home/yuzy/synthetic/grtm/example/crust1.0_45/10.grn.{i}")[0]

    d = data2[:, idxs[i]]
    t = np.arange(len(d)) * 0.02
    ax.plot(t, d, alpha=0.5, c="r")
    #ax.plot(t, tr.data, alpha=0.5, c="b")

    #print(np.max(np.abs(d)), np.max(np.abs(tr.data)))
    #print(np.argmax(np.abs(tr.data)))

axs = [fig.add_subplot((gs[i, 1])) for i in range(6)]
for i in range(10):
    ax = axs[i%6]

    d = data2[:, i]
    t = np.arange(len(d)) * 0.02
    ax.plot(t, d, alpha=0.5, c="r")
    #print(np.max(np.abs(d)), np.max(np.abs(tr.data)))
    #print(np.argmax(np.abs(tr.data)))
#plt.xlim((5, 15))
plt.savefig("img4.pdf")