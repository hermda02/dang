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

ad = [1.6444e-2, 8.39923e-3, 1.19689e-3, 5.890824e-2, 2.9665593e-1]

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
plt.close()


A_d_1 = np.loadtxt('sim_data_020_A_d.dat')                                                
like = np.loadtxt('sim_data_020_likelihood_A_d.dat')                                    
plt.plot(A_d_1,like)
plt.axvline(A_d_1[50],color='r', label='Sampled')
plt.axvline(ad[0],color='g',label='Nominal')
plt.title(r'Likelihood slice: $A_d$ at 20 GHz',size=20)
plt.xlabel(r'$A_d$',size=15)
plt.ylabel('Likelihood',size=15)
plt.savefig('likelihood_slice_A_d_020.png',dpi=300,bbox_inches='tight') 
plt.close()

A_d_2 = np.loadtxt('sim_data_030_A_d.dat')                                                
like = np.loadtxt('sim_data_030_likelihood_A_d.dat')                                    
plt.plot(A_d_2,like)
plt.axvline(A_d_2[50],color='r', label='Sampled')
plt.axvline(ad[1],color='g',label='Nominal')
plt.title(r'Likelihood slice: $A_d$ at 30 GHz',size=20)
plt.xlabel(r'$A_d$',size=15)
plt.ylabel('Likelihood',size=15)
plt.savefig('likelihood_slice_A_d_030.png',dpi=300,bbox_inches='tight') 
plt.close()

A_d_3 = np.loadtxt('sim_data_045_A_d.dat')                                                
like = np.loadtxt('sim_data_045_likelihood_A_d.dat')                                    
plt.plot(A_d_3,like)
plt.axvline(A_d_3[50],color='r', label='Sampled')
plt.axvline(ad[2],color='g',label='Nominal')
plt.title(r'Likelihood slice: $A_d$ at 45 GHz',size=20)
plt.xlabel(r'$A_d$',size=15)
plt.ylabel('Likelihood',size=15)
plt.savefig('likelihood_slice_A_d_045.png',dpi=300,bbox_inches='tight') 
plt.close()

A_d_4 = np.loadtxt('sim_data_070_A_d.dat')                                                
like = np.loadtxt('sim_data_070_likelihood_A_d.dat')                                    
plt.plot(A_d_4,like)
plt.axvline(A_d_4[50],color='r', label='Sampled')
plt.axvline(ad[3],color='g',label='Nominal')
plt.title(r'Likelihood slice: $A_d$ at 70 GHz',size=20)
plt.xlabel(r'$A_d$',size=15)
plt.ylabel('Likelihood',size=15)
plt.savefig('likelihood_slice_A_d_070.png',dpi=300,bbox_inches='tight') 
plt.close()

A_d_5 = np.loadtxt('sim_data_100_A_d.dat')                                                
like  = np.loadtxt('sim_data_100_likelihood_A_d.dat')                                    
plt.plot(A_d_5,like)
plt.axvline(A_d_5[50],color='r', label='Sampled')
plt.axvline(ad[4],color='g',label='Nominal')
plt.title(r'Likelihood slice: $A_d$ at 100 GHz',size=20)
plt.xlabel(r'$A_d$',size=15)
plt.ylabel('Likelihood',size=15)
plt.savefig('likelihood_slice_A_d_100.png',dpi=300,bbox_inches='tight') 
plt.close()
