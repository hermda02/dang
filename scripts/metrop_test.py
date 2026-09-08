import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

mpl.rcParams['text.usetex'] = True

missval = -1.6375e30

files = ['accept.dat', 'prob.dat', 'chi.dat', 'temps.dat']

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

USAGE = f"Usage: python3 {sys.argv[0]} [--help] | -accept -beta -chisq -prob"

def plot() -> None:
    command = sys.argv[1:]
    if not command:
        raise SystemExit(USAGE)

    for i in command:
        if i == '--help':
            raise SystemExit(USAGE)
        elif i == '-accept':
            plot_accept()
        elif i == '-beta':
            trace_beta()
        elif i == '-chisq':
            trace_chisq()
        elif i == '-prob':
            plot_prob()
        else:
            raise SystemExit(USAGE)


if __name__ == "__main__":
    plot()
