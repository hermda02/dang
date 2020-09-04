import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import stats

# mpl.rcParams['text.usetex'] = True

plt.rc('text', usetex=True)

dust_Q = np.loadtxt('dust_Q_amplitudes.dat')
dust_U = np.loadtxt('dust_U_amplitudes.dat')

freq   = [28.4, 44.1, 70.3, 33.0, 40.6]
labels = ['bp\_030','bp\_044','bp\_070','WMAP\_K','WMAP\_Q']

print(f"Data files have {np.shape(dust_Q)[0]} + iterations.")
iter = input("Which iteration would you like to plot? ")

x = np.linspace(0.5*np.min(freq),1.5*np.max(freq),100)

def power_law(amp,ref,x):
    return amp*(x/ref)**(-3.1)

fig, ax = plt.subplots(figsize=(6,6))
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

ax.plot(x,power_law(dust_Q[num,3],freq[3],x),color='k',label='synch')
ax.set_title('Dust Q Amplitudes',size=20)
ax.set_ylim([0.5e-3,1e-1])
ax.set_xlim([0.9*np.min(freq),1.1*np.max(freq)])
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
