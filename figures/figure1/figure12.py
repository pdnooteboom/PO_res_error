import numpy as np
from matplotlib.pylab import plt
import seaborn as sns

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

import ot
import ot.plot

from mpl_toolkits import mplot3d

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

sns.set(style='whitegrid', font_scale=1.2)

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

#%%
from matplotlib import cm
import matplotlib.colors as colors

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

fig = plt.figure(figsize=(10,8))

#ax = plt.subplot(211)
ax = plt.subplot2grid((9,3), (0,0), colspan=3, rowspan=5, projection='3d')
ax.invert_zaxis()

x = np.linspace(-40, 25, 10)
y = np.linspace(-40, 25, 10)
X, Y = np.meshgrid(x, y)
Z = 1.5+np.random.rand(X.shape[0],X.shape[1])/3.
#ax.contour3D(X, Y, Z, color='k', zorder=-10)
ax.plot_surface(X, Y, Z, cmap=truncate_colormap(cm.Reds, 0.3, 1),
                       linewidth=0, antialiased=False, vmin=1.5, vmax=1.75, alpha=0.3, zorder=-100)
#ax.contourf(X, Y, Z, cmap=cm.coolwarm,
#                        antialiased=False,vmin=0.5, vmax=1.8, alpha=0.7, zorder=-100)



plt.xlim(-40,25)
plt.ylim(-40,25)
ax.set_xlabel('$^{\circ}$E', fontsize=18)
ax.set_ylabel('$^{\circ}$N', fontsize=18)
ax.set_zlabel('depth (km)', fontsize=18)


leng = 20
xs = np.linspace(xL+20, x1[1], leng) + 10*np.sin(np.linspace(0,4*np.pi,leng))
ys = np.linspace(yL, y1[1], leng) + 1*np.cos(np.linspace(0,4*np.pi,leng))
zs = np.linspace(0.9, 0, leng)+ 0.1*np.sin(np.linspace(0,2*np.pi,leng))

ax.plot(xs, ys, zs,':', color='k', linewidth = 2, zorder = 10)

#a = Arrow3D([xL+20, x1[0]], [yL+4, y1[0]], 
#            [1, 0], mutation_scale=20, 
#            lw=3, arrowstyle="->", color="k", zorder=10)
#ax.add_artist(a)

ax.scatter3D(x1,y1, color=firstcloudcolor,alpha=1, s=50, label='first distribution')
ax.scatter3D(xL+17,yL, [0.9], color='k', marker = 'P', s=255, label = 'release location', zorder=10)

ax.zaxis.set_ticks([1,0])
ax.zaxis.set_ticklabels([1,0])


#ax.tick_params(axis='x',labelbottom=False, labelleft=False, colors='red', width=0) 
#ax.tick_params(axis='y',labelbottom=False, labelleft=False, colors='red', width=0) 
ax.set_yticks([])
ax.set_xticks([])
                   

plt.title('(a)', fontsize=18)


#%%

ax = plt.subplot2grid((9,3), (6,2), rowspan=3)
plt.plot(np.concatenate((x1[np.newaxis,:],x2[np.where(G0>0)[1]][np.newaxis,:]), axis=0), 
         np.concatenate((y1[np.newaxis,:],y2[np.where(G0>0)[1]][np.newaxis,:]), axis=0),
         '--', color='k', alpha=0.8)
plt.scatter(x1,y1, color=firstcloudcolor,alpha=1, label='first distribution')
plt.scatter(x2,y2, color=secondcloudcolor,alpha=1, marker = '+', s=60, label='second distribution')
plt.scatter(xL,yL, color='k', marker = 'P', s=135, label = 'release location')


plt.xlabel('$^{\circ}$E', fontsize=18)
ax.tick_params(labelbottom=False, labelleft=False) 
plt.title('(d)', fontsize=18)
plt.xlim(-40,40)
plt.ylim(-40,40)

plt.legend(bbox_to_anchor=(2.10, 1.5))


#%%
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

#ax = plt.subplot(338)
ax = plt.subplot2grid((9,3), (6,1), rowspan=3)
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
plt.title('(c)', fontsize=18)
plt.xlim(-40,40)
plt.ylim(-40,40)

#%%
#ax = plt.subplot(337)
ax = plt.subplot2grid((9,3), (6,0), rowspan=3)
plt.plot(np.concatenate((np.full(n,xL)[np.newaxis,:], x1[np.newaxis,:]),axis=0), np.concatenate((np.full(n,yL)[np.newaxis,:], y1[np.newaxis,:]),axis=0), color='red', zorder=1)
plt.scatter(x1,y1, color=firstcloudcolor,alpha=1, label='first cloud',zorder=2)#,color='firebrick'
plt.scatter(xL,yL, color='k', marker = 'P', s=135, label = 'release location')
plt.xlabel('$^{\circ}$E', fontsize=18)
plt.ylabel('$^{\circ}$N', fontsize=18)
ax.tick_params(labelbottom=False, labelleft=False) 
plt.title('(b)', fontsize=18)
plt.xlim(-40,40)
plt.ylim(-40,40)



plt.savefig('figure1.pdf', bbox_inches="tight")
plt.show()