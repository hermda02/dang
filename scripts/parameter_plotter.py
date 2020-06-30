import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import corner
from scipy import stats
import sys

# mpl.rcParams['text.usetex'] = True

plt.rc('text', usetex=True)

A_d_mean = 11.4378598319838
A_s_mean = 1.69279291160855
b_s_mean = -3.10000
means    = [A_s_mean, A_d_mean, b_s_mean]

# For beautiful maps
planck  = np.loadtxt('/home/daniel/graduate_school/masters/tools/Planck_color.txt')/255.
pl_cmap = colors.LinearSegmentedColormap.from_list('planck',planck)

files = ['pixel_100_A_s.dat', 'pixel_100_A_d.dat', 'pixel_100_beta_s.dat']

data = []
labels = []
ranges = np.empty((3,2))

for file in files:
    data.append(np.loadtxt(file))
    try:
        labels.append((file.replace("pixel_100_","")).replace(".dat",""))
    except:
        try:
            labels.append(file.replace(".dat",""))
        except:
            labels.append(file)
            
labels = [label.replace("_", "\_") for label in labels]
            
chi = np.loadtxt('total_chisq.dat')

iterations = len(data[0])

for i in range(3):
    ranges[i] = (np.mean(data[i])-2.0*np.std(data[i]),np.mean(data[i])+2.0*np.std(data[i]))

a_s = data[0]
a_d = data[1]
b_s = data[2]
x   = np.linspace(1,iterations,iterations)

def getminmax(x):
    mx = np.percentile(x, 97.5)
    mn = np.percentile(x, 2.5)  
    return (mn, mx)

def hjornet(burnin):

    samples = np.vstack(data).T
    samples = samples[burnin:]
    ndim    = np.shape(samples)[1]
    num     = np.shape(samples)[0]
    figure  = corner.corner(samples, labels=labels, bins=int(np.sqrt(num)), range=ranges)
    axes    = np.array(figure.axes).reshape((ndim, ndim))
    col     = ['r','g','b']

    for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(means[i],color=col[i])

    for xi in range(ndim):
        for yi in range(xi):
            ax = axes[xi, yi]
            ax.axvline(means[yi],color=col[yi])
            ax.axhline(means[xi],color=col[xi])

    plt.savefig('corner_plot.png',dpi=300,bbox_inches='tight')

def a_s_plot():
    binnum = int(np.sqrt(iterations))
    
    fig,axes = plt.subplots(1,2,figsize=(8,4))
    fig.tight_layout(pad=2.0)
    axes[0].plot(x,a_s)
    axes[0].set_title(r'$A_s$ Trace',size=15)
    axes[0].set_xlabel('Gibbs Iteration',size=10)
    axes[0].set_ylabel(r'$A_s$',size=15)
    axes[0].axhline(y=A_s_mean,c='k')
    axes[1].hist(a_s,bins=binnum)
    axes[1].set_xlabel(r'$A_s$',size=10)
    axes[1].set_ylabel('Count',size=15)
    axes[1].axvline(x=A_s_mean,c='k')
    plt.savefig('a_s_trace_hist',dpi=150,bbox_inches='tight')
    #plt.show()
    plt.close()

def b_s_plot():
    binnum = int(np.sqrt(iterations))
    
    fig,axes = plt.subplots(1,2,figsize=(8,4))
    fig.tight_layout(pad=2.0)
    axes[0].plot(x,b_s)
    axes[0].set_title(r'$\beta_s$ Trace',size=15)
    axes[0].set_xlabel('Gibbs Iteration',size=10)
    axes[0].set_ylabel(r'$\beta_s$',size=15)
    axes[0].axhline(y=b_s_mean,c='k')
    axes[1].hist(b_s,bins=binnum)
    axes[1].set_xlabel(r'$\beta_s$',size=10)
    axes[1].set_ylabel('Count',size=15)
    axes[1].axvline(x=b_s_mean,c='k')
    plt.savefig('beta_s_trace_hist',dpi=150,bbox_inches='tight')
    #plt.show()
    plt.close()

def chisq_plot():
    binnum = int(np.sqrt(iterations))
    
    fig,axes = plt.subplots(1,2,figsize=(8,4))
    fig.tight_layout(pad=2.0)
    axes[0].plot(x,chi)
    axes[0].set_yscale('log')
    axes[0].set_title('Trace',size=15)
    axes[0].set_ylabel(r'$\chi^2$',size=15)
    axes[0].set_xlabel('Gibbs Iteration',size=10)
    axes[1].hist(chi,bins=binnum)
    axes[1].set_ylabel('Count',size=15)
    axes[1].set_xlabel(r'$\chi^2$',size=10)
    plt.savefig('chisquare',dpi=150,bbox_inches='tight')
    plt.close()
    # plt.show()

def a_b_s_histo2d():
    
    a_s_dif = (ranges[0][1] - ranges[0][0])/(2*np.sqrt(iterations))
    b_s_dif = (ranges[2][1] - ranges[2][0])/(2*np.sqrt(iterations))

    xedges = np.arange(ranges[0][0], ranges[0][1], a_s_dif)
    yedges = np.arange(ranges[2][0], ranges[2][1], b_s_dif)

    d, xedges, yedges = np.histogram2d(a_s, b_s, bins=(xedges,yedges)) 
    d1 = d.T

    x1, y1 = np.meshgrid(xedges,yedges)

    fig = plt.figure(figsize=(8,6))

    plt.xlim(np.min(a_s),np.max(a_s))
    plt.ylim(np.min(b_s),np.max(b_s))
    plt.pcolormesh(x1, y1, d1,vmin=0,vmax=5,cmap=plt.get_cmap('hot_r'))
    plt.axvline(x=A_s_mean,color='green')
    plt.axhline(y=-3.1,color='k')
    plt.xlabel('Dust Amplitude',size=15)
    plt.ylabel(r'$\beta_s$',size=15)
    im1 = plt.pcolormesh(x1, y1, d1,vmin=0,vmax=5,cmap=plt.get_cmap('hot_r'))
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar_ax.tick_params(labelsize=20)
    cbar = fig.colorbar(im1,cax=cbar_ax)
    cbar.set_label('Pixel Count',size=18)
    plt.savefig('beta_s_A_s_2dhist',dpi=150,bbox_inches='tight')
    #plt.show()
    plt.close()


USAGE = f"Usage: python3 {sys.argv[0]} [--help] | -chisq -A_s -beta -corner -A_beta"

def plot() -> None:
    command = sys.argv[1:]
    if not command:
        raise SystemExit(USAGE)

    for i in command:
        if (i == '--help'):
            SystemExit(USAGE)
        if (i == '-chisq'):
            chisq_plot()
        if (i == '-A_s'):
            a_s_plot()
        if (i == '-beta'):
            b_s_plot()
        if (i == '-A_beta'):
            a_b_s_histo2d()
        if (i == '-corner'):
            burn = int(input(f"Burn-in value? (Total of {iterations} iterations): "))
            hjornet(burn)
        else:
            SystemExit(USAGE)

plot()
