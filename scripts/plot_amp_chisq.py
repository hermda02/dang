import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import stats

# mpl.rcParams['text.usetex'] = True

plt.rc('text', usetex=True)

# For beautiful maps
planck  = np.loadtxt('/home/daniel/graduate_school/masters_thesis/tools/Planck_color.txt')/255.
pl_cmap = colors.LinearSegmentedColormap.from_list('planck',planck)

a_s = np.loadtxt('pixel_150_A_s.dat')
b_s = np.loadtxt('pixel_150_beta_s.dat')
chi = np.loadtxt('total_chisq.dat')

iterations = len(a_s)

x   = np.linspace(1,iterations,iterations)

# plt.plot(x,a_s)
# plt.show()

# plt.plot(x,b_s)
# plt.show()

# plt.plot(x,chi)
# plt.yscale('log')
# plt.show()

xedges = np.arange(-3, -1.99, 0.01)
yedges = np.arange(-3.3, -2.795, 0.005)

d, xedges, yedges = np.histogram2d(a_s, b_s, bins=(xedges,yedges)) 
d1 = d.T
x1, y1 = np.meshgrid(xedges,yedges)

fig = plt.figure(figsize=(8,6))

plt.pcolormesh(x1, y1, d1, cmap=plt.get_cmap('hot_r'))
plt.xlim(-3.0,-2.0)
plt.ylim(-3.3,-2.8)
im1 = plt.pcolormesh(x1, y1, d1,vmin=0,vmax=500,cmap=plt.get_cmap('hot_r'))
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar_ax.tick_params(labelsize=20)
cbar = fig.colorbar(im1,cax=cbar_ax)
cbar.set_label('Pixel Count',size=18)
plt.show()