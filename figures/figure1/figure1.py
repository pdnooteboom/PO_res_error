import numpy as np
from matplotlib.pylab import plt
import seaborn as sns

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

import ot
import ot.plot

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

sns.set(style='darkgrid', font_scale=1.2)

surfacecolor = 'dodgerblue'
firstcloudcolor = 'k'
secondcloudcolor = 'forestgreen'

#%%
xL = -30; yL = -30;
sigma = 9
sigma2 = 8
bias = 10

res = 3

con = 3
con2 = 32

n = 8

np.random.seed(1)

x1 = np.random.normal(xL+bias,sigma2,n) + 12*con
x2 = np.random.normal(xL,sigma,n)+14 

y1 = np.random.normal(yL,sigma2+2,n) + 16
y2 = np.random.normal(yL+bias,sigma,n)+con2 

#Define OT
M = ot.dist(np.concatenate((x1[:,np.newaxis],y1[:,np.newaxis]), axis=1), np.concatenate((x2[:,np.newaxis],y2[:,np.newaxis]), axis=1))
M /= M.max()

G0 = ot.emd(np.ones((n,)) / n, np.ones((n,)) / n, M)



sns.set_style("dark")


fig = plt.figure(figsize=(10,2))
ax = fig.add_subplot(133)
plt.plot(np.concatenate((x1[np.newaxis,:],x2[np.where(G0>0)[1]][np.newaxis,:]), axis=0), 
         np.concatenate((y1[np.newaxis,:],y2[np.where(G0>0)[1]][np.newaxis,:]), axis=0),
         '--', color='k', alpha=0.8)
plt.scatter(x1,y1, color=firstcloudcolor,alpha=1, label='first cloud')
plt.scatter(x2,y2, color=secondcloudcolor,alpha=1, marker = '+', s=60, label='second cloud')
plt.scatter(xL,yL, color='k', marker = 'P', s=135, label = 'release location')


plt.xlabel('$^{\circ}$E', fontsize=18)
ax.tick_params(labelbottom=False, labelleft=False) 
plt.title('(c)', fontsize=18)
plt.xlim(-40,40)
plt.ylim(-40,40)

plt.legend(bbox_to_anchor=(1.92, 1.05))

step = 8
xs, ys = np.mgrid[-44:48:step,-44:48:step]
vs = np.ones(xs.shape,dtype=bool)

xss = xs[:,0]
yss = ys[0]

boxes = []
for i in range(len(x1)):
    nex = find_nearest_index(xss,x1[i])
    ney = find_nearest_index(yss,y1[i])
    if(xss[nex]<x1[i]):
        lbx = xss[nex]
        hbx = xss[nex+1]
    else:
        lbx = xss[nex-1]
        hbx = xss[nex]
    if(yss[ney]<y1[i]):
        lby = yss[ney]
        hby = yss[ney+1]
    else:
        lby = yss[ney-1]
        hby = yss[ney]

    boxes.append(Rectangle((lbx, lby), hbx-lbx, hby-lby))
    
    ix = int(x1[i]/step)
    iy = int(y1[i]/step)
    vs[ix, iy] = 0
    
gridboxes = []
for i in range(xs.shape[0]):
    for j in range(xs.shape[1]):
        gridboxes.append(Rectangle((xs[i,j], ys[i,j]), step, step))

zs = np.ma.array(~vs, mask=vs).astype(int)*1000

ax = fig.add_subplot(132)
pc2 = PatchCollection(gridboxes, facecolor='lightgray', alpha=0.6,
                         edgecolor='lightgray')
ax.add_collection(pc2)
pc = PatchCollection(boxes, facecolor=surfacecolor, alpha=1,
                         edgecolor=surfacecolor)
ax.add_collection(pc)
plt.scatter(x1,y1, color=firstcloudcolor,alpha=1, label='first cloud')#,color='firebrick'
plt.scatter(xL,yL, color='k', marker = 'P', s=135, label = 'release location')
plt.xlabel('$^{\circ}$E', fontsize=18)
ax.tick_params(labelbottom=False, labelleft=False) 
plt.title('(b)', fontsize=18)
plt.xlim(-40,40)
plt.ylim(-40,40)

ax = fig.add_subplot(131)
plt.plot(np.concatenate((np.full(n,xL)[np.newaxis,:], x1[np.newaxis,:]),axis=0), np.concatenate((np.full(n,yL)[np.newaxis,:], y1[np.newaxis,:]),axis=0), color='red', zorder=1)
plt.scatter(x1,y1, color=firstcloudcolor,alpha=1, label='first cloud',zorder=2)#,color='firebrick'
plt.scatter(xL,yL, color='k', marker = 'P', s=135, label = 'release location')
plt.xlabel('$^{\circ}$E', fontsize=18)
plt.ylabel('$^{\circ}$N', fontsize=18)
ax.tick_params(labelbottom=False, labelleft=False) 
plt.title('(a)', fontsize=18)
plt.xlim(-40,40)
plt.ylim(-40,40)



plt.savefig('figure1.pdf', bbox_inches="tight")
plt.show()