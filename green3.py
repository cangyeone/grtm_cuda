import ctypes 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as grid 
import numpy as np 
import obspy 
r_earth = 6371 #地球半径
def readtxt(path):
    f = open(path, "r") 
    data = []
    for line in f.readlines():
        data.append([float(r) for r in line.strip().split(",")[:-1] if len(r)>0])
    return np.array(data)
def readpar(path, flatten=False):
    f = open(path, "r") 
    data = []
    r = r_earth 
    tdata = []
    for line in f.readlines():
        #0.0000,6.1000,3.5500,2.7400,1.000e+03,5.000e+02
        sline = [float(r) for r in line.strip().split(",") if len(r)>0]
        tdata.append(sline) 
    tdata2 = []
    for idx, line in enumerate(tdata):
        h, vp, vs, rho, qp, qs = line
        if idx == len(tdata)-1:
            sh = 0 
        else:
            sh = tdata[idx+1][0]-tdata[idx][0]
        r -= sh 
        
        if flatten:
            s = r_earth / (r + 0.5 * sh)
        else:
            s = 1 
        sh *= s 
        vs *= s 
        vp *= s 
        tdata2.append([sh, vp, vs, rho, qp, qs])
        data.append([sh, vp, vs, rho, qp, qs])
    h = 0 
    for idx, line in enumerate(tdata2):
        if idx == 0:
            h = 0 
        else:
            h += tdata2[idx-1][0]
        data[idx][0] = h   
    return data 
f_iname = "demo/input.txt"
f_oname_w = "demo/oname_w.txt"
f_oname_t = "demo/oname_t.txt"
model = "demo/crust1.0.txt"
flatten = True 
data = readpar(model, flatten)
edep = 25 #震源深度 
rdep = 0.0  #接收点深度
nlayer = len(data) 
dists = 1000 #震中距
t0 = 56.19 # 
delta = 0.5 
log2_npts = 10  
src_type = 0 
n_ptam = 7 
taper = 0.3 

if flatten:
    edep = r_earth * np.log(r_earth / (r_earth - edep))
    rdep = r_earth * np.log(r_earth / (r_earth - rdep))


ofile = open(f_iname, "w")
strs = f"""{edep} 
{rdep}
{nlayer}
"""
ofile.write(strs)
for h, vp, vs, rho, qp, qs in data:
    ofile.write(f"{h},{vp},{vs},{rho},{qp},{qs}\n")

tags = f"""{delta}
{dists} 
1
1000
33.33 
{log2_npts}
{t0}
{taper}
{n_ptam}
{src_type}

"""
ofile.write(tags)
ofile.close()

func = ctypes.cdll.LoadLibrary("build/libpopm.so")
iname = ctypes.c_char_p(f_iname.encode()) 
oname_w = ctypes.c_char_p(f_oname_w.encode())
oname_t = ctypes.c_char_p(f_oname_t.encode())
func.green(iname, oname_w, oname_t)

data = readtxt(f_oname_t)
fig = plt.figure(1, figsize=(18, 18), dpi=300)
gs = grid.GridSpec(6, 1)


idxs = [4, 1, 5, 3, 0, 2]

for i in range(6):
    ax = fig.add_subplot(gs[i, 0])
    tr = obspy.read(f"../grtm_fortran/example/crust1.0_{int(edep)}/{dists}.grn.{i}")[0]
    d = data[:, idxs[i]]
    if i == 2:d = d * 0.0 
    t = np.arange(len(d)) * delta 
    #print(tr.stats)
    d2 = tr.data 
    t2 = np.arange(len(d2)) * tr.stats.delta 
    ax.plot(t, d, alpha=0.5, c="r")
    ax.plot(t2, d2, alpha=0.5, c="b")
plt.savefig("demo/demo3.png")
    

EALL = 0
for i in range(6):
    tr = obspy.read(f"../grtm_fortran/example/crust1.0_{int(edep)}/{dists}.grn.{i}")[0]
    d = data[:, idxs[i]]
    if i == 2:d = d * 0.0 
    d2 = tr.data 
    ad1 = np.abs(d) 
    ad2 = np.abs(d2) 
    if np.sum(ad2)==0:continue 
    E = np.sum(np.abs(ad1-ad2))/(np.sum(ad2)+np.sum(ad1))
    EALL += E 
print("ERROR", EALL/5) 