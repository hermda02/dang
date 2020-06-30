import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.stats import norm
from scipy.optimize import curve_fit
from decimal import Decimal

mpl.rcParams['text.usetex'] = True

missval = -1.6375e30

# For beautiful maps
planck  = np.loadtxt('/home/daniel/graduate_school/masters/tools/Planck_color.txt')/255.
pl_cmap = colors.LinearSegmentedColormap.from_list('planck',planck)

files  = ['accept.dat'. 'prob.dat', 'chi.dat', 'temps.dat']

#-----------------------------------------

def plot_accept():

    accept = np.loadtxt(files[0])

    niter = np.linspace(1,len(accept),len(accept))

    plt.plot(niter,accept,color='orange')
    plt.title('Acceptance Rate',size=20)
    plt.xlabel('Iteration Number',size=20)
    plt.ylabel(r'$N_{accept}/N_{\rm iter}$',size=20)
    plt.savefig('accept_trace',dpi=300,bbox_inches='tight')
    plt.close()

def plot_prob():

    prob = np.loadtxt(files[1])

    niter = np.linspace(1,len(prob),len(prob))

    plt.plot(niter,prob,color='orange')
    plt.title('Acceptance Probability',size=20)
    plt.xlabel('Iteration Number',size=20)
    plt.ylabel('Probability',size=20)
    plt.savefig('accept_prob',dpi=300,bbox_inches='tight')
    plt.close()
    
def trace_beta():

    beta = np.loadtxt(files[3])

    niter = np.linspace(1,len(beta),len(beta))

    plt.plot(niter,beta,color='orange')
    plt.title(r'Trace of $\beta_s$',size=20)
    plt.xlabel('Iteration Number',size=20)
    plt.ylabel(r'$\beta_s$',size=20)
    plt.savefig('beta_trace',dpi=300,bbox_inches='tight')
    plt.close()

def trace_chisq():

    chisq = np.loadtxt(files[2])

    niter = np.linspace(1,len(chisq),len(chisq))

    plt.plot(niter,chisq,color='orange')
    plt.title(r'Trace of $\chi^2$',size=20)
    plt.xlabel('Iteration Number',size=20)
    plt.yscale('log')
    plt.ylabel(r'$\chi^2$',size=20)
    plt.savefig('chisq_trace',dpi=300,bbox_inches='tight')
    plt.close()

USAGE = f"Usage: python3 {sys.argv[0]} [--help] |  -accept -beta -chisq -prob"

def plot() -> None:
    commander = sys.argv[1:]
    if not command:
        raise SystemExit(USAGE)

    for i in command:
        if (i == '--help'):
            SystemExit(USAGE)
        if (i == '-accept'):
            plot_accept()
        if (i == '-beta'):
            trace_beta()
        if (i == '-chisq'):
            trace_chisq()
        if (i == '-prob'):
            plot_prob()
        else:
            SystemExit(USAGE)
