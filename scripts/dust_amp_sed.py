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

# print()
direct = input("Which directory do you want to plot from? ")

dust_Q = np.loadtxt(direct+'dust_Q_amplitudes.dat')
dust_U = np.loadtxt(direct+'dust_U_amplitudes.dat')

freq = [28.4, 44.1, 70.1, 70.1, 70.1, 22.8, 33.0, 40.6, 40.6, 60.8, 60.8]#, 143.0]   

print(f"Data files have {np.shape(dust_Q)[0]} + iterations.")
iter = input("Which iteration would you like to plot? ")
num = int(iter)-1

plt.scatter(freq,dust_Q[num,:])
plt.title('Dust Q Amplitudes',size=20)
plt.ylim([0,0.15])
plt.xlabel('Frequency [GHz]',size=15)
plt.ylabel('Dust Amplitude',size=15)
plt.savefig(direct+'amplitudes_dust_Q_'+iter+'.png',dpi=300,bbox_inches='tight')
plt.close()

plt.scatter(freq,dust_U[num,:])
plt.title('Dust U Amplitudes',size=20)
plt.ylim([0,0.15])
plt.xlabel('Frequency [GHz]',size=15)
plt.ylabel('Dust Amplitude',size=15)
plt.savefig(direct+'amplitudes_dust_U_'+iter+'.png',dpi=300,bbox_inches='tight')
plt.close()
