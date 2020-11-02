import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import stats

# mpl.rcParams['text.usetex'] = True

plt.rc('text', usetex=True)

ame_spec = np.loadtxt('../../ame_spectra/spdust2_cnm.dat')

ame_spec = ame_spec.T

# Used to shift the spdust spectrum                                                                                                                                                                                
def sdust(nu, Asd, nu_p, nuref=22.):
    nu_ref = nuref*1e9
    nu_p0 = 30.*1e9

    fnu, f = np.loadtxt("../../ame_spectra/spdust2_cnm.dat", unpack=True)
    fnu *= 1e9
    # MAKE SURE THAT THESE ARE BOTH IN 1e9                                                                                                                                                                         
    scale = nu_p0/nu_p

    f = np.interp(scale*nu*1e9, fnu, f)
    f0 = np.interp(scale*nu_ref, scale*nu*1e9, f) # Value of s at nu_0                                                                                                                                             
    s_sd = Asd*f/f0
    return s_sd

def nearest(array,val):
    array = np.asarray(array)
    idx = (np.abs(array-val)).argmin()
    return idx

ame_spec2 = sdust(ame_spec[0],1.0,17.0e9)

dust_Q = np.loadtxt('dust_Q_amplitudes.dat')
dust_U = np.loadtxt('dust_U_amplitudes.dat')

freq   = [28.4, 44.1, 70.3, 22.8, 40.6]
labels = ['bp\_030','bp\_044','bp\_070','WMAP\_K','WMAP\_Q']

print(f"Data files have {np.shape(dust_Q)[0]} + iterations.")
iter = input("Which iteration would you like to plot? ")

x = np.linspace(0.5*np.min(freq),1.5*np.max(freq),100)

def power_law(amp,ref,x):
    return amp*(x/ref)**(-3.1)


fig, ax = plt.subplots(figsize=(8,6))
if type(iter) == "all":
    for j in range(len(freq)):
        if dust_Q[0][j] == 0.0:
            continue
        else:
            for i in range(len(dust_Q[0])):
                ax.plot(freq[j],dust_Q[i,j])
else:
    num = int(iter)-1
    for j in range(len(freq)):
        if dust_Q[num][j] == 0.0:
            continue
        else:
            ax.scatter(freq[j],dust_Q[num,j],label=labels[j])


z = nearest(ame_spec[0],22.8)

#print(x)
#print(ame_spec[0][x])
#print(ame_spec2[x])
#exit()
ame_spec2 = ame_spec2*(dust_Q[num,3]/ame_spec2[z])
ax.plot(x,power_law(dust_Q[num,3],freq[3],x),color='k',label='synch')
#ax.plot(ame_spec[0],ame_spec2,'k--',label='spdust2 cnm')
ax.set_title('Dust Q Amplitudes',size=20)
ax.set_ylim([0.5e-3,0.5e-1])
ax.set_xlim([20,60])
ax.set_xlabel('Frequency [GHz]',size=15)
ax.set_yscale('log')
ax.set_ylabel('Dust Amplitude',size=15)
ax.legend()
plt.savefig('amplitudes_dust_Q_'+iter+'.png',dpi=300,bbox_inches='tight')
plt.close()

#plt.scatter(freq,dust_U[num,:])
#plt.title('Dust U Amplitudes',size=20)
#plt.ylim([0,0.15])
#plt.xlabel('Frequency [GHz]',size=15)
#plt.ylabel('Dust Amplitude',size=15)
#plt.savefig(direct+'amplitudes_dust_U_'+iter+'.png',dpi=300,bbox_inches='tight')
#plt.close()
