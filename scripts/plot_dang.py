import numpy as np
import healpy as hp
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import stats
import sys
import os
import re

hp.disable_warnings()

def gaussian(x,amp,mu,std):
    return amp*np.exp(-0.5*((x-mu)/std)**2)

def read_params(filename):
    global mu, sd
    labels = []
    freq   = []
    with open(filename,'r') as infile:
        for line in infile:
            if line.startswith('NUMBAND'):
                numbands = int(line.split('=')[1])
            if line.startswith('NUMGIBBS'):
                numgibbs = int(line.split('=')[1][:5])
            if line.startswith('MASKFILE'):
                maskfile = str.strip(line.split('=')[1])

    with open(filename,'r') as infile:
        for line in infile:
            if line.startswith('COMP_PRIOR_GAUSS_BETA'):
                if "MEAN" in line:
                    mu = float(str.strip(line.split('=')[1]))
                if "STD" in line:
                    sd = float(str.strip(line.split('=')[1]))
    blabs = []
    bfreq = []
    for band in range(numbands):
        blabs.append('BAND_LABEL'+str(band+1).zfill(3))
        bfreq.append('BAND_FREQ'+str(band+1).zfill(3))        
    for band in range(numbands):
         with open(filename,'r') as infile:
            for line in infile:
                if line.startswith(blabs[band]):
                    name = str.strip(line.split('=')[1])
                    labels.append(name)
                if line.startswith(bfreq[band]):
                    fre  = str.strip(line.split('=')[1])
                    freq.append(float(fre))
    return labels, freq, numgibbs, maskfile

def _init_():
    plt.rc('text', usetex=True)

    global dir, names, freq, num_samp, labels, num_bands, iterations, files
    global maskfile

    try:    
        files = os.listdir()   
        thing = [file for file in files if 'param' in file]
        names, freq, num_samp, maskfile = read_params(thing[0])
        labels     = [name.replace("'","") for name in names]
        num_bands  = len(freq)
        iterations = num_samp
        maskfile   = '../'+maskfile.replace("'","")
        print(mu,sd)
        print(labels)
        print(freq)
    except Exception as e:
        print(e)
        print("Input which directory you wish to point to (../ automatically included)")
        exit()

def load_data():
    global diag, parameters, chiQ, chiU, chisq
    global param_Q, param_U
    global asQ, adQ, bsQ, adU, asU, bsU, x
    global ranges_Q, ranges_U
    global data_Q, data_U
    global mask

    global amplitudes
    diag = []
    parameters = []

    amplitudes = np.loadtxt('dust_353_Q_amplitudes.dat')

    chiQ = np.loadtxt('total_chisq_Q.dat')    
    chiU = np.loadtxt('total_chisq_U.dat') 

    chisq = np.sqrt(chiQ**2 + chiU**2)

    mask = hp.read_map(maskfile)
    
    for file in files:
        if file.startswith('pixel') and file.endswith('dat'):
            diag.append(file)
            pixel = re.search('_(.+?)_',file)[1]

    diag_Q = [dia for dia in diag if "Q" in dia]
    diag_U = [dia for dia in diag if "U" in dia]

    # for target in diag:
    #     parameters.append((target.replace("pixel_"+pixel+"_","")).replace(".dat",""))


    parameters = [par.replace("pixel_"+pixel+"_","").replace(".dat","") for par in diag]
    parameters = [par.replace("_","\_") for par in parameters]
    param_Q = [par for par in parameters if "Q" in par]
    param_U = [par for par in parameters if "U" in par]
    param_Q.sort()
    param_U.sort()

    data_Q = []
    data_U = []

    ranges_Q = np.empty((3,2))
    ranges_U = np.empty((3,2))


    for file in diag_Q:
        data_Q.append(np.loadtxt(file))
        if "beta" in file:
            bsQ = np.loadtxt(file)
        if "A_s" in file:
            asQ = np.loadtxt(file)
        if "A_d" in file:
            adQ = np.loadtxt(file)
    for file in diag_U:
        data_U.append(np.loadtxt(file))
        if "beta" in file:
            bsU = np.loadtxt(file)
        if "A_s" in file:
            asU = np.loadtxt(file)
        if "A_d" in file:
            adU = np.loadtxt(file)

    for i in range(3):
        ranges_Q[i] = (np.mean(data_Q[i])-4.5*np.std(data_Q[i]),np.mean(data_Q[i])+4.5*np.std(data_Q[i]))
        ranges_U[i] = (np.mean(data_U[i])-4.5*np.std(data_U[i]),np.mean(data_U[i])+4.5*np.std(data_U[i]))

    # bsQ = data_Q[0]
    # asQ = data_Q[1]
    # adQ = data_Q[2]
    # adU = data_U[0]
    # bsU = data_U[1]
    # asU = data_U[2]
    x   = np.linspace(1,iterations,iterations)
    
def correlate_dust_amps(burnin):
    a_ame = np.loadtxt('dust_353_Q_amplitudes.dat')

    label = [lab.replace("_", "\_") for lab in labels]

    a_ame  = a_ame[burnin:][:]

    stats_list = label
    df = pd.DataFrame(a_ame, columns=stats_list)
    for i in range(np.shape(a_ame)[1]):
        if a_ame[0][i] == 0.0:
            del df[stats_list[i]]

    df.insert(5,r'$\chi^2_{_{Q}}$',chiQ,True)
    df.insert(6,r'$\chi^2_{_{U}}$',chiU,True)

    corr = df.corr()

    f, ax = plt.subplots(figsize=(10, 9))
    plt.title(r'$a_{\rm AME}$ Correlation')
    sns.heatmap(corr, mask=np.zeros_like(corr, dtype=np.bool),annot=True, 
    cmap=sns.diverging_palette(220, 10, as_cmap=True), square=True,  
            ax=ax, vmin=-1.0, vmax=1.0)

    plt.savefig('a_ame_corr_plot.pdf', bbox_inches='tight')
    # plt.show()

def trace_all(pol):

    sim = False#True
    
    if sim:
        sim_synch_030 = hp.read_map('../data/test_data/new_sims/synch/synch_030_n0064_rj.fits',field=(1,2))
        noise_sim_030 = hp.read_map('../data/test_data/new_sims/noise/sim_030_noise_sim_nominal_n0064.fits',field=(1,2))
        true_val_Q  = sim_synch_030[0][23000]
        true_val_U  = sim_synch_030[1][23000]
        noise_val_Q = noise_sim_030[0][23000]
        noise_val_U = noise_sim_030[1][23000]

        true_chi_q  = np.loadtxt('../sim0.50_v64/total_chisq_Q.dat')
        true_chi_u  = np.loadtxt('../sim0.50_v64/total_chisq_U.dat')
        
    binnum = int(np.sqrt(iterations))

    betax = np.linspace(-3.5,-2.5,100)
    
    if pol == "Q":

        # A_s
        fig,axes = plt.subplots(4,2,figsize=(8,16))
        fig.tight_layout(pad=2.0)
        axes[0][0].plot(x,asQ)
        axes[0][0].set_xlabel('Gibbs Iteration',size=10)
        axes[0][0].set_ylabel(r'$A_s$',size=15)
        axes[0][1].hist(asQ,bins=binnum)
        axes[0][1].set_xlabel(r'$A_s$',size=10)
        axes[0][1].set_ylabel('Count',size=15)

        # A_d
        axes[1][0].plot(x,adQ)
        axes[1][0].set_ylabel(r'$A_d$',size=15)
        axes[1][1].hist(adQ,bins=binnum)
        axes[1][1].set_xlabel(r'$A_d$',size=10)
        axes[1][1].set_ylabel('Count',size=15)

        # Beta_s
        axes[2][0].plot(x,bsQ)
        axes[2][0].set_ylabel(r'$\beta_s$',size=15)
        axes[2][1].hist(bsQ,bins=binnum)
        axes[2][1].plot(betax,gaussian(betax,40,-3.1,0.1))
        axes[2][1].set_xlim([np.min(bsQ),np.max(bsQ)])
        axes[2][1].set_xlabel(r'$\beta_s$',size=10)
        axes[2][1].set_ylabel('Count',size=15)

        # Chisq
        axes[3][0].plot(x,chiQ)
        axes[3][0].set_ylabel(r'$\chi^2$',size=15)
        axes[3][0].set_xlabel('Gibbs Iteration',size=10)
        axes[3][1].hist(chiQ,bins=binnum)
        axes[3][1].set_ylabel('Count',size=15)
        axes[3][1].set_xlabel(r'$\chi^2$',size=10)
        plt.savefig('all_trace_Q',dpi=150,bbox_inches='tight')
        plt.close()

    if pol == "U":
        binnum = int(np.sqrt(iterations))

        # A_s
        fig,axes = plt.subplots(4,2,figsize=(8,16))
        fig.tight_layout(pad=2.0)
        axes[0][0].plot(x,asU)
        axes[0][0].set_xlabel('Gibbs Iteration',size=10)
        axes[0][0].set_ylabel(r'$A_s$',size=15)
        axes[0][1].hist(asU,bins=binnum)
        axes[0][1].set_xlabel(r'$A_s$',size=10)
        axes[0][1].set_ylabel('Count',size=15)

        # A_d
        axes[1][0].plot(x,adU)
        axes[1][0].set_ylabel(r'$A_d$',size=15)
        axes[1][1].hist(adU,bins=binnum)
        axes[1][1].set_xlabel(r'$A_d$',size=10)
        axes[1][1].set_ylabel('Count',size=15)

        # Beta_s
        axes[2][0].plot(x,bsU)
        axes[2][0].set_ylabel(r'$\beta_s$',size=15)
        axes[2][1].hist(bsU,bins=binnum)
        axes[2][1].set_xlabel(r'$\beta_s$',size=10)
        axes[2][1].set_ylabel('Count',size=15)

        # Chisq
        axes[3][0].plot(x,chiU)
        axes[3][0].set_ylabel(r'$\chi^2$',size=15)
        axes[3][0].set_xlabel('Gibbs Iteration',size=10)
        axes[3][1].hist(chiU,bins=binnum)
        axes[3][1].set_ylabel('Count',size=15)
        axes[3][1].set_xlabel(r'$\chi^2$',size=10)
        plt.savefig('all_trace_U',dpi=150,bbox_inches='tight')
        plt.close()

def a_d_trace():
    binnum = int(np.sqrt(iterations))

    band = int(input(f"Which band? {labels} "))
    amps = np.loadtxt('dust_353_Q_amplitudes.dat')

    print(np.shape(amps))

    samples = np.shape(amps)[0]

    x_arr = np.arange(samples)
    
    data_ad = amps.T[band]
    
    print(np.mean(data_ad),np.std(data_ad))
    
    fig,axes = plt.subplots(1,2,figsize=(8,4))
    fig.tight_layout(pad=2.0)
    axes[0].plot(x_arr,data_ad)
    axes[0].set_ylabel(r'$A_d$',size=15)
    axes[0].set_xlabel('Gibbs Iteration',size=10)
    axes[1].hist(data_ad,bins=binnum)
    axes[1].set_xlabel(r'$A_d$',size=10)
    axes[1].set_ylabel('Count',size=15)
    plt.savefig('Ad_trace_'+labels[band],dpi=150,bbox_inches='tight')
    plt.close()
    # if pol == "U":
    #     fig,axes = plt.subplots(1,2,figsize=(8,4))
    #     fig.tight_layout(pad=2.0)
    #     axes[0].plot(x,adU)
    #     axes[0].set_xlabel('Gibbs Iteration',size=10)
    #     axes[0].set_ylabel(r'$A_d$',size=15)
    #     axes[1].hist(adQ,bins=binnum)
    #     axes[1].set_xlabel(r'$A_d$',size=10)
    #     axes[1].set_ylabel('Count',size=15)
    #     plt.savefig('Ad_trace_Q',dpi=150,bbox_inches='tight')
    #     plt.close()

def a_s_trace(pol):

    if pol == "Q":
        binnum = int(np.sqrt(iterations))
    
        fig,axes = plt.subplots(1,2,figsize=(8,4))
        fig.tight_layout(pad=2.0)
        axes[0].plot(x,asQ)
        axes[0].set_title(r'$A_s$ Q Trace',size=15)
        axes[0].set_xlabel('Gibbs Iteration',size=10)
        axes[0].set_ylabel(r'$A_s$',size=15)
        # axes[0].axhline(y=A_s_mean,c='k')
        axes[1].hist(asQ,bins=binnum)
        axes[1].set_xlabel(r'$A_s$',size=10)
        axes[1].set_ylabel('Count',size=15)
        # axes[1].axvline(x=A_s_mean,c='k')
        plt.savefig('a_s_trace_hist_Q',dpi=150,bbox_inches='tight')
        #plt.show()
        plt.close()

    if pol == "U":
        binnum = int(np.sqrt(iterations))
    
        fig,axes = plt.subplots(1,2,figsize=(8,4))
        fig.tight_layout(pad=2.0)
        axes[0].plot(x,asU)
        axes[0].set_title(r'$A_s$ U Trace',size=15)
        axes[0].set_xlabel('Gibbs Iteration',size=10)
        axes[0].set_ylabel(r'$A_s$',size=15)
        # axes[0].axhline(y=A_s_mean,c='k')
        axes[1].hist(asU,bins=binnum)
        axes[1].set_xlabel(r'$A_s$',size=10)
        axes[1].set_ylabel('Count',size=15)
        # axes[1].axvline(x=A_s_mean,c='k')
        plt.savefig('a_s_trace_hist_U',dpi=150,bbox_inches='tight')
        #plt.show()
        plt.close()

def b_s_trace(pol):
    y = np.linspace(mu-3*sd,mu+3*sd,1000)
    if pol == "Q":
        binnum = int(np.sqrt(iterations))
    
        fig,axes = plt.subplots(1,2,figsize=(8,4))
        fig.tight_layout(pad=2.0)
        axes[0].plot(x,bsQ)
        axes[0].set_title(r'$\beta_s$ Trace',size=15)
        axes[0].set_xlabel('Gibbs Iteration',size=10)
        axes[0].set_ylabel(r'$\beta_s$',size=15)
        # axes[0].axhline(y=b_s_mean,c='k')
        axes[1].hist(bsQ,bins=binnum,density=True)
        axes[1].set_xlabel(r'$\beta_s$',size=10)
        axes[1].set_ylabel('Count',size=15)
        axes[1].plot(y,stats.norm.pdf(y,mu,sd))
        # axes[1].axvline(x=b_s_mean,c='k')
        plt.savefig('beta_s_trace_hist_Q',dpi=150,bbox_inches='tight')
        #plt.show()
        plt.close()

    if pol == "U":
        binnum = int(np.sqrt(iterations))
    
        fig,axes = plt.subplots(1,2,figsize=(8,4))
        fig.tight_layout(pad=2.0)
        axes[0].plot(x,bsU)
        axes[0].set_title(r'$\beta_s$ Trace',size=15)
        axes[0].set_xlabel('Gibbs Iteration',size=10)
        axes[0].set_ylabel(r'$\beta_s$',size=15)
        # axes[0].axhline(y=b_s_mean,c='k')
        axes[1].hist(bsU,bins=binnum)
        axes[1].set_xlabel(r'$\beta_s$',size=10)
        axes[1].set_ylabel('Count',size=15)
        axes[1].plot(y,stats.norm.pdf(y,mu,sd))
        # axes[1].axvline(x=b_s_mean,c='k')
        plt.savefig('beta_s_trace_hist_Q',dpi=150,bbox_inches='tight')
        #plt.show()
        plt.close()
        
        
def chisq_trace(pol):
    if pol == "Q":
        binnum = int(np.sqrt(iterations))
    
        fig,axes = plt.subplots(1,2,figsize=(8,4))
        fig.tight_layout(pad=2.0)
        axes[0].plot(x,chiQ)
        axes[0].set_yscale('log')
        axes[0].set_title('Trace',size=15)
        axes[0].set_ylabel(r'$\chi^2$',size=15)
        axes[0].set_xlabel('Gibbs Iteration',size=10)
        axes[1].hist(chiQ,bins=binnum)
        axes[1].set_ylabel('Count',size=15)
        axes[1].set_xlabel(r'$\chi^2$',size=10)
        plt.savefig('chisquare_Q',dpi=150,bbox_inches='tight')
        plt.close()
        # plt.show()

    if pol == "U":
        binnum = int(np.sqrt(iterations))
    
        fig,axes = plt.subplots(1,2,figsize=(8,4))
        fig.tight_layout(pad=2.0)
        axes[0].plot(x,chiU)
        axes[0].set_yscale('log')
        axes[0].set_title('Trace',size=15)
        axes[0].set_ylabel(r'$\chi^2$',size=15)
        axes[0].set_xlabel('Gibbs Iteration',size=10)
        axes[1].hist(chiU,bins=binnum)
        axes[1].set_ylabel('Count',size=15)
        axes[1].set_xlabel(r'$\chi^2$',size=10)
        plt.savefig('chisquare_U',dpi=150,bbox_inches='tight')
        plt.close()
        # plt.show()

def hjornet(burnin,pol):
    import corner
    
    if pol == "Q":
        samples = np.vstack(data_Q).T
        samples = samples[burnin:]
        ndim    = np.shape(samples)[1]
        num     = np.shape(samples)[0]
        val     = np.mean(samples,axis=0)
        std     = np.std(samples,axis=0)
        figure  = corner.corner(samples, labels=param_Q,range=ranges_Q) # bins=20, 
        axes    = np.array(figure.axes).reshape((ndim, ndim))
        col     = ['r','g','b']
        
        for i in range(ndim):
            ax = axes[i,i]
            ax.axvline(val[i],color=col[i])
            ax.axvline(val[i]+std[i],color=col[i],linestyle='--')
            ax.axvline(val[i]-std[i],color=col[i],linestyle='--')

        for yi in range(ndim):
            for xi in range(yi):
                ax = axes[yi,xi]
                ax.axvline(val[xi],color=col[xi])
                ax.axhline(val[yi],color=col[yi])
            
        # plt.show()
        plt.savefig('corner_plot_Q.png',dpi=300,bbox_inches='tight')

    if pol == "U":
        samples = np.vstack(data_Q).T
        samples = samples[burnin:]
        ndim    = np.shape(samples)[1]
        num     = np.shape(samples)[0]
        val     = np.mean(samples,axis=0)
        std     = np.std(samples,axis=0)
        figure  = corner.corner(samples, labels=parameters,range=ranges_Q) # bins=20, 
        axes    = np.array(figure.axes).reshape((ndim, ndim))
        col     = ['r','g','b']
        
        for i in range(ndim):
            ax = axes[i,i]
            ax.axvline(val[i],color=col[i])
            ax.axvline(val[i]+std[i],color=col[i],linestyle='--')
            ax.axvline(val[i]-std[i],color=col[i],linestyle='--')

        for yi in range(ndim):
            for xi in range(yi):
                ax = axes[yi,xi]
                ax.axvline(val[xi],color=col[xi])
                ax.axhline(val[yi],color=col[yi])
            
        plt.savefig('corner_plot_U.png',dpi=300,bbox_inches='tight')

def beta_chisq(pol):
    y = np.linspace(mu-3*sd,mu+3*sd,1000)
    sim = True
    
    if sim:
        sim_synch_030 = hp.read_map('../data/test_data/new_sims/synch/synch_030_n0064_rj.fits',field=(1,2))
        noise_sim_030 = hp.read_map('../data/test_data/new_sims/noise/sim_030_noise_sim_nominal_n0064.fits',field=(1,2))
        true_val_Q  = sim_synch_030[0][23000]
        true_val_U  = sim_synch_030[1][23000]
        noise_val_Q = noise_sim_030[0][23000]
        noise_val_U = noise_sim_030[1][23000]

        true_chi_q  = np.loadtxt('../sim0.50_v64/total_chisq_Q.dat')
        true_chi_u  = np.loadtxt('../sim0.50_v64/total_chisq_U.dat')
        
    binnum = int(np.sqrt(iterations))

    if pol == "Q":

        # A_s
        fig,axes = plt.subplots(2,2,figsize=(8,8))
        fig.tight_layout(pad=2.0)

        axes[0][0].plot(x,bsQ)
        axes[0][0].set_ylabel(r'$\beta_s$',size=15)
        axes[0][0].axhline(-3.1,linestyle='--',color='k',label='True beta')
        axes[0][1].hist(bsQ,bins=binnum)
        axes[0][1].set_xlabel(r'$\beta_s$',size=10)
        axes[0][1].set_ylabel('Count',size=15)
        axes[0][1].axvline(-3.1,linestyle='--',color='k',label='True beta')

        # Chisq
        axes[1][0].plot(x,chiQ)
        axes[1][0].axhline(true_chi_q,color='k',linestyle='--',label='true solution')
        axes[1][0].legend()
        axes[1][0].set_ylabel(r'$\chi^2$',size=15)
        axes[1][0].set_xlabel('Gibbs Iteration',size=10)
        axes[1][1].hist(chiQ,bins=binnum)
        axes[1][1].set_ylabel('Count',size=15)
        axes[1][1].set_xlabel(r'$\chi^2$',size=10)
        plt.savefig('beta_chi_trace_Q',dpi=150,bbox_inches='tight')
        plt.close()

        
        # chi_dif_Q = (1.10*np.max(chiQ) - np.min(chiQ))/(2*np.sqrt(iterations))
        # b_s_dif_Q = (ranges_Q[2][1] - ranges_Q[2][0])/(2*np.sqrt(iterations))

        # print(chi_dif_Q)

        # xedges = np.arange(np.min(chiQ), np.max(chiQ),chi_dif_Q)
        # yedges = np.arange(ranges_Q[2][0], ranges_Q[2][1], b_s_dif_Q)

        # print(xedges)

        # d, xedges, yedges = np.histogram2d(chiQ, bsQ, bins=(xedges,yedges)) 
        # d1 = d.T

        # x1, y1 = np.meshgrid(xedges,yedges)


        # fig, axes = plt.subplots(1,2,figsize=(12,6),sharey=True)

        # axes[0].plot(chiQ,bsQ)
        # axes[0].set_xlim(np.min(chiQ),np.max(chiQ))
        # axes[0].set_ylim(np.min(bsQ),np.max(bsQ))
        # axes[0].set_ylabel(r'$\beta_s$',size=15)
        # axes[0].set_xlabel(r'$\chi^2$',size=15)

        # axes[1].set_xlim(np.min(chiQ),np.max(chiQ))
        # axes[1].set_ylim(np.min(bsQ),np.max(bsQ))
        # axes[1].pcolormesh(x1, y1, d1,vmin=0,vmax=binnum,cmap=plt.get_cmap('hot_r'))

        # axes[1].set_xlabel(r'$\chi^2$',size=15)
        # # axes[1].set_ylabel(r'$\beta_s$ Q',size=15)
        # im1 = plt.pcolormesh(x1, y1, d1,vmin=0,vmax=binnum,cmap=plt.get_cmap('hot_r'))
        # cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        # cbar_ax.tick_params(labelsize=20)
        # cbar = fig.colorbar(im1,cax=cbar_ax)
        # cbar.set_label('Pixel Count',size=18)

        # plt.savefig('beta_chisq_Q',dpi=150,bbox_inches='tight')
        # plt.close()

    if pol == "U":

        # A_s
        fig,axes = plt.subplots(2,2,figsize=(8,8))
        fig.tight_layout(pad=2.0)

        axes[0][0].plot(x,bsU)
        axes[0][0].set_ylabel(r'$\beta_s$',size=15)
        axes[0][0].axhline(-3.1,linestyle='--',color='k',label='True beta')
        axes[0][1].hist(bsU,bins=binnum)
        axes[0][1].set_xlabel(r'$\beta_s$',size=10)
        axes[0][1].set_ylabel('Count',size=15)
        axes[0][1].axvline(-3.1,linestyle='--',color='k',label='True beta')

        # Chisq
        axes[1][0].plot(x,chiU)
        axes[1][0].axhline(true_chi_u,color='k',linestyle='--',label='true solution')
        axes[1][0].legend()
        axes[1][0].set_ylabel(r'$\chi^2$',size=15)
        axes[1][0].set_xlabel('Gibbs Iteration',size=10)
        axes[1][1].hist(chiU,bins=binnum)
        axes[1][1].set_ylabel('Count',size=15)
        axes[1][1].set_xlabel(r'$\chi^2$',size=10)
        plt.savefig('beta_chi_trace_U',dpi=150,bbox_inches='tight')
        plt.close()

        
    # if pol == "U":

    #     fig, axes = plt.subplots(figsize=(6,6))

    #     axes.plot(chiU,bsU)
    #     axes.set_ylabel(r'$\beta_s$',size=15)
    #     axes.set_xlabel(r'$\chi^2$',size=15)
    #     plt.savefig('beta_chisq_U',dpi=150,bbox_inches='tight')
    #     plt.close()
        
def a_b_s_histo2d(pol):
    if pol == "Q":
     
        a_s_dif_Q = (ranges_Q[0][1] - ranges_Q[0][0])/(2*np.sqrt(iterations))
        b_s_dif_Q = (ranges_Q[2][1] - ranges_Q[2][0])/(2*np.sqrt(iterations))

        xedges = np.arange(ranges_Q[0][0], ranges_Q[0][1], a_s_dif_Q)
        yedges = np.arange(ranges_Q[2][0], ranges_Q[2][1], b_s_dif_Q)

        d, xedges, yedges = np.histogram2d(asQ, bsQ, bins=(xedges,yedges)) 
        d1 = d.T

        x1, y1 = np.meshgrid(xedges,yedges)

        fig = plt.figure(figsize=(8,6))

        plt.xlim(np.min(asQ),np.max(asQ))
        plt.ylim(np.min(bsQ),np.max(bsQ))
        plt.pcolormesh(x1, y1, d1,vmin=0,vmax=5,cmap=plt.get_cmap('hot_r'))
        # plt.axvline(x=A_s_mean,color='green')
        # plt.axhline(y=-3.1,color='k')
        plt.xlabel('Synch Amplitude Q',size=15)
        plt.ylabel(r'$\beta_s$ Q',size=15)
        im1 = plt.pcolormesh(x1, y1, d1,vmin=0,vmax=5,cmap=plt.get_cmap('hot_r'))
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar_ax.tick_params(labelsize=20)
        cbar = fig.colorbar(im1,cax=cbar_ax)
        cbar.set_label('Pixel Count',size=18)
        plt.savefig('beta_s_A_s_2dhist_Q',dpi=150,bbox_inches='tight')
        #plt.show()
        plt.close()

    if pol == "U":
     
        a_s_dif_U = (ranges_U[0][1] - ranges_U[0][0])/(2*np.sqrt(iterations))
        b_s_dif_U = (ranges_U[2][1] - ranges_U[2][0])/(2*np.sqrt(iterations))

        xedges = np.arange(ranges_U[0][0], ranges_U[0][1], a_s_dif_U)
        yedges = np.arange(ranges_U[2][0], ranges_U[2][1], b_s_dif_U)

        d, xedges, yedges = np.histogram2d(asU, bsU, bins=(xedges,yedges)) 
        d1 = d.T

        x1, y1 = np.meshgrid(xedges,yedges)

        fig = plt.figure(figsize=(8,6))

        plt.xlim(np.min(asU),np.max(asU))
        plt.ylim(np.min(bsU),np.max(bsU))
        plt.pcolormesh(x1, y1, d1,vmin=0,vmax=5,cmap=plt.get_cmap('hot_r'))
        # plt.axvline(x=A_s_mean,color='green')
        # plt.axhline(y=-3.1,color='k')
        plt.xlabel('Dust Amplitude U',size=15)
        plt.ylabel(r'$\beta_s$ U',size=15)
        im1 = plt.pcolormesh(x1, y1, d1,vmin=0,vmax=5,cmap=plt.get_cmap('hot_r'))
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar_ax.tick_params(labelsize=20)
        cbar = fig.colorbar(im1,cax=cbar_ax)
        cbar.set_label('Pixel Count',size=18)
        plt.savefig('beta_s_A_s_2dhist_U',dpi=150,bbox_inches='tight')
        #plt.show()
        plt.close()

def getminmax(x):
    mx = np.percentile(x, 97.5)
    mn = np.percentile(x, 2.5)  
    return (mn, mx)

def return_mean_map(list,outfile):
    print('Creating '+str(outfile))
    nside = 64
    npix  = nside*nside*12
    samples = len(list)

    print(samples)

    maps = np.empty((npix,samples))
    out_map = np.empty(npix)
    
    maps = maps.T
    
    for i in range(samples):
        maps[i] = hp.read_map(list[i],verbose=False)

    maps = maps.T

    for i in range(npix):
        out_map[i] = np.mean(maps[i][:])

    hp.write_map(str(outfile),out_map)

def return_std_map(list,outfile):
    print('Creating '+str(outfile))
    nside = 64
    npix  = nside*nside*12
    samples = len(list)

    print(samples)

    maps = np.empty((npix,samples))
    out_map = np.empty(npix)
    
    maps = maps.T
    
    for i in range(samples):
        maps[i] = hp.read_map(list[i],verbose=False)

    maps = maps.T

    for i in range(npix):
        out_map[i] = np.std(maps[i][:])

    hp.write_map(outfile,out_map)

def make_mean_maps():

    res_Q = []
    res_U = []

    for j in range(len(files)):
        file = files[j]
        for i in range(len(labels)):
            if file.startswith(labels[i]+'_residual_Q'):
                if file.endswith('.fits'):
                    res_Q.append(file)
    for j in range(len(files)):
        file = files[j]
        for i in range(len(labels)):
            if file.startswith(labels[i]+'_residual_Q'):
                if file.endswith('.fits'):
                    res_U.append(file)

    for i in range(num_bands):
        maps = []
        for res in res_Q:   
            if labels[i] in res:
                maps.append(res)
        return_mean_map(maps,labels[i]+'_residual_Q_mean.fits')
        maps = []
        for res in res_U:   
            if labels[i] in res:
                maps.append(res)
        return_mean_map(maps,labels[i]+'_residual_U_mean.fits')

    components = ['chisq_Q', 'chisq_U', 'bp_030_synch_amplitude_Q', 'bp_030_synch_amplitude_U', 'synch_beta_Q']
    comps = ['chisq_Q', 'chisq_U', 'synch_Q_030', 'synch_U_030', 'synch_beta']

    for i in range(len(components)):
        maps = []
        for file in files:
            if file.startswith(components[i]) and file.endswith('.fits'):
                maps.append(file)
        return_mean_map(maps,comps[i]+'_mean.fits')        
        if comps[i] == 'synch_Q':
            return_std_map(maps,comps[i]+'_bp30_std.fits')
        if comps[i] == 'synch_U':
            return_std_map(maps,comps[i]+'_bp30_std.fits')
        if comps[i] == 'synch_beta':
            return_std_map(maps,comps[i]+'_std.fits')

def make_synch_diff_maps(mask):

    joint_synch_Q = hp.read_map('synch_mean.fits',field=1)
    joint_synch_U = hp.read_map('synch_mean.fits',field=2)

    beta = hp.read_map('synch_beta_mean.fits',field=1)

    npix = len(beta)

    joint_Q_30 = np.empty(npix)
    joint_U_30 = np.empty(npix)
    # joint_Q_spass = np.empty(npix)
    # joint_U_spass = np.empty(npix)
    
    for i in range(npix):
        if mask[i] == 0:
            joint_Q_30[i] = hp.UNSEEN#(30./44.1)**beta[i] * joint_synch_Q[i]
            joint_U_30[i] = hp.UNSEEN#(30./44.1)**beta[i] * joint_synch_U[i]
            # joint_Q_spass[i] = hp.UNSEEN#(2.305/44.1)**beta[i] * joint_synch_Q[i]
            # joint_U_spass[i] = hp.UNSEEN#(2.305/44.1)**beta[i] * joint_synch_U[i]
        else:
            joint_Q_30[i] = (30./28.4)**beta[i] * joint_synch_Q[i]
            joint_U_30[i] = (30./28.4)**beta[i] * joint_synch_U[i]
            # joint_Q_spass[i] = (2.305/44.1)**beta[i] * joint_synch_Q[i]
            # joint_U_spass[i] = (2.305/44.1)**beta[i] * joint_synch_U[i]
    
    bp_030_Q      = hp.read_map('../data/BP_synch_Q_n0064.fits')
    bp_030_U      = hp.read_map('../data/BP_synch_U_n0064.fits')

    npipe_30_Q    = hp.read_map('../data/npipe6v20_comm_synch_n0064_60arc_Q_rc1.fits')
    npipe_30_U    = hp.read_map('../data/npipe6v20_comm_synch_n0064_60arc_U_rc1.fits')

    # spass_Q       = hp.read_map('../data/spass_rmrot_n0064_ring_masked.fits',field=1)
    # spass_U       = hp.read_map('../data/spass_rmrot_n0064_ring_masked.fits',field=2)

    bp_min_joint_Q    = bp_030_Q - joint_Q_30
    bp_min_joint_U    = bp_030_U - joint_U_30

    np_min_joint_Q    = npipe_30_Q - joint_Q_30
    np_min_joint_U    = npipe_30_U - joint_U_30

    # spass_min_joint_Q = spass_Q - joint_Q_spass
    # spass_min_joint_U = spass_U - joint_U_spass

    
    for i in range(npix):
        if mask[i] == 0:
            bp_min_joint_Q[i] = hp.UNSEEN#    = bp_030_Q - joint_Q_30
            bp_min_joint_U[i] = hp.UNSEEN#    = bp_030_U - joint_U_30
            np_min_joint_Q[i] = hp.UNSEEN#    = npipe_30_Q - joint_Q_30
            np_min_joint_U[i] = hp.UNSEEN#    = npipe_30_U - joint_U_30
            # spass_min_joint_Q[i] = hp.UNSEEN# = spass_Q - joint_Q_spass
            # spass_min_joint_U[i] = hp.UNSEEN# = spass_U - joint_U_spass

    hp.write_map('npipe_minus_joint_Q_60arcmin_n0064.fits',np_min_joint_Q)
    hp.write_map('npipe_minus_joint_U_60arcmin_n0064.fits',np_min_joint_U)
    hp.write_map('BP_minus_joint_Q_60arcmin_n0064.fits',bp_min_joint_Q)
    hp.write_map('BP_minus_joint_U_60arcmin_n0064.fits',bp_min_joint_U)
    # hp.write_map('spass_minus_joint_Q_60arcmin_n0064.fits',spass_min_joint_Q)
    # hp.write_map('spass_minus_joint_U_60arcmin_n0064.fits',spass_min_joint_U)

def trace_mean(data):

    print(np.shape(data))
    n  = np.shape(data.T)[1]
    print(n)

    print(np.linspace(0,num_bands-1,num_bands))
    print(freq)
    # print(np.vstack(data[:5,:]))
    band = int(input(""))

    x  = np.linspace(1,n,n)
    xm = np.linspace(2,n,n-1)


    # fig, ax = plt.subplots(num_bands,1)

    # for band in range(num_bands):
    #     means = np.zeros(n-1)
    #     stds  = np.zeros(n-1)
    #     for i in range(1,n-1):
    #         means[i] = np.mean(data[:i,band])
    #         stds[i]  = np.std(data[:i,band])

    #     ax[band].plot(x,data[:,band])
    #     ax[band].errorbar(xm,means,yerr=stds,color='k')
    means = np.zeros(n-1)
    stds  = np.zeros(n-1)
    for i in range(1,n-1):
        means[i] = np.mean(data[:i,band])
        stds[i]  = np.std(data[:i,band])

    plt.plot(x,data[:,band])
    plt.errorbar(xm,means,yerr=stds,color='k')
    plt.show()

USAGE = f"Usage: python3 {sys.argv[0]} [paramfile] [option]\n Option list: \n -Ad \n -As \n -As_beta \n -beta \n -beta_chi \n -chisq \n -corner \n -correlate \n -mean_maps \n -synch_diff \n -trace_all"

def plot() -> None:
    command = sys.argv[1:]
    if not command:
        raise SystemExit(USAGE)
    _init_()
    load_data()
    for i in command:
        if (i == '--help'):
            SystemExit(USAGE)
        elif (i == '-chisq'):
            pols = str(input("Q or U? "))
            chisq_trace(pols)
        elif (i == '-As'):
            pols = str(input("Q or U? "))
            a_s_trace(pols)
        elif (i == '-beta'):
            pols = str(input("Q or U? "))
            b_s_trace(pols)
        elif (i == '-Ad'):
            a_d_trace()
        elif (i == '-trace_all'):
            pols = str(input("Q or U? "))
            trace_all(pols)
        elif (i == '-beta_chi'):
            pols = str(input("Q or U? "))
            beta_chisq(pols)
        elif (i == '-As_beta'):
            pols = str(input("Q or U? "))
            a_b_s_histo2d(pols)
        elif (i == '-corner'):
            burn = int(input(f"Burn-in value? (Total of {iterations} iterations): "))
            pols = str(input("Q or U? "))
            hjornet(burn,pols)
        elif (i == '-correlate'):
            burn = int(input(f"Burn-in value? (Total of {iterations} iterations): "))
            correlate_dust_amps(burn)
        elif (i == '-mean_maps'):
            make_mean_maps()
        elif (i == '-synch_diff'):
            make_synch_diff_maps(mask)
        elif (i == '-trace_mean'):
            par = 'dust'#input(f"Which parameter? \n {parameters}")
            # if par == 'dust':
                # print(amplitudes)
            trace_mean(amplitudes)
        else:
            SystemExit(USAGE)

plot()
