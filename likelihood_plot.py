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

A_s = np.loadtxt('A_s.dat')                                                
like = np.loadtxt('likelihood_A_s.dat')                                    
plt.plot(A_s,like)
plt.axvline(A_s[50],color='r',label='Sampled')
plt.axvline(23.1543446,color='g',label='Nominal')
plt.title(r'Likelihood slice: $A_s$ at pixel 350',size=20)
plt.xlabel(r'$A_s$',size=15)
plt.ylabel('Likelihood',size=15)
plt.legend(loc='best')
plt.savefig('likelihood_slice_A_s.png',dpi=300,bbox_inches='tight') 
plt.show()
