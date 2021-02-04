import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os
import sys

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
        maps[i] = hp.read_map('../'+dir+'/'+list[i],verbose=False)

    maps = maps.T

    for i in range(npix):
        out_map[i] = np.mean(maps[i][:])

    hp.write_map('../'+dir+'/'+str(outfile),out_map)

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
        maps[i] = hp.read_map('../'+dir+'/'+list[i],verbose=False)

    maps = maps.T

    for i in range(npix):
        out_map[i] = np.std(maps[i][:])

    hp.write_map('../'+dir+'/'+outfile,out_map)
    
def read_params(filename):
    labels = []
    freq   = []
    with open(filename,'r') as infile:
        for line in infile:
            if line.startswith('NUMBAND'):
                numbands = int(line.split('=')[1])
            if line.startswith('NUMGIBBS'):
                numgibbs = int(line.split('=')[1][:5])
                
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
    return labels, freq, numgibbs

try:    
    dir = sys.argv[1] #str(input('Which directory to mean? '))
    files = os.listdir('../'+dir)
    names, freq, num_samp = read_params('../'+dir+'/param_'+dir+'.txt')
    labels = [name.replace("'","") for name in names]
    num_bands = len(freq)
    print(labels)
    print(freq)

except:
    print("Input which directory you wish to point to (../ automatically included)")
    exit()

res_1_Q = []
res_2_Q = []
res_3_Q = []
res_4_Q = []
res_5_Q = []
res_6_Q = []
res_7_Q = []
res_8_Q = []

res_1_U = []
res_2_U = []
res_3_U = []
res_4_U = []
res_5_U = []
res_6_U = []
res_7_U = []
res_8_U = []

synch_Qs = []
synch_Us = []
chi_Qs   = []
chi_Us   = []

beta_ss  = []

for file in files:
    
    # Residuals Q
    
    if file.startswith(labels[0]+'_residual_Q'):
        if file.endswith('.fits'):
            res_1_Q.append(file)

    if file.startswith(labels[1]+'_residual_Q'):
        if file.endswith('.fits'):
            res_2_Q.append(file)
    
    if file.startswith(labels[2]+'_residual_Q'):
        if file.endswith('.fits'):
            res_3_Q.append(file)

    if file.startswith(labels[3]+'_residual_Q'):
        if file.endswith('.fits'):
            res_4_Q.append(file)

    if file.startswith(labels[4]+'_residual_Q'):
        if file.endswith('.fits'):
            res_5_Q.append(file)

    if file.startswith(labels[5]+'_residual_Q'):
        if file.endswith('.fits'):
            res_6_Q.append(file)

    if file.startswith(labels[6]+'_residual_Q'):
        if file.endswith('.fits'):
            res_7_Q.append(file)

#    if file.startswith(labels[7]+'_residual_Q'):
#        if file.endswith('.fits'):
#            res_8_Q.append(file)

    # Residuals U
    
    if file.startswith(labels[0]+'_residual_U'):
        if file.endswith('.fits'):
            res_1_U.append(file)

    if file.startswith(labels[1]+'_residual_U'):
        if file.endswith('.fits'):
            res_2_U.append(file)
    
    if file.startswith(labels[2]+'_residual_U'):
        if file.endswith('.fits'):
            res_3_U.append(file)

    if file.startswith(labels[3]+'_residual_U'):
        if file.endswith('.fits'):
            res_4_U.append(file)

    if file.startswith(labels[4]+'_residual_U'):
        if file.endswith('.fits'):
            res_5_U.append(file)
            
    if file.startswith(labels[5]+'_residual_U'):
        if file.endswith('.fits'):
            res_6_U.append(file)

    if file.startswith(labels[6]+'_residual_U'):
        if file.endswith('.fits'):
            res_7_U.append(file)
            
    # if file.startswith(labels[7]+'_residual_U'):
    #     if file.endswith('.fits'):
    #         res_8_U.append(file)
            
    # Synch maps
            
    if file.startswith('bp_044_synch_amplitude_Q'):
        if file.endswith('.fits'):
            synch_Qs.append(file)
            
    if file.startswith('bp_044_synch_amplitude_U'):
        if file.endswith('.fits'):
            synch_Us.append(file)

    # Synch beta
    if file.startswith('synch_beta_Q'):
        if file.endswith('.fits'):
            beta_ss.append(file)
            
    # Chisquare maps

    if file.startswith('chisq_Q'):
        if file.endswith('.fits'):
            chi_Qs.append(file)

    if file.startswith('chisq_U'):
        if file.endswith('.fits'):
            chi_Us.append(file)


return_mean_map(chi_Qs,'chisq_Q_mean.fits')
return_mean_map(chi_Us,'chisq_U_mean.fits')
return_mean_map(synch_Qs,'synch_Q_mean.fits')
return_mean_map(synch_Us,'synch_U_mean.fits')

return_mean_map(beta_ss,'synch_beta_mean.fits')

return_std_map(synch_Qs,'synch_Q_std.fits')
return_std_map(synch_Us,'synch_U_std.fits')

return_std_map(beta_ss,'synch_beta_std.fits')

return_mean_map(res_1_Q,labels[0]+'_residual_Q_mean.fits')
return_mean_map(res_2_Q,labels[1]+'_residual_Q_mean.fits')
return_mean_map(res_3_Q,labels[2]+'_residual_Q_mean.fits')
return_mean_map(res_4_Q,labels[3]+'_residual_Q_mean.fits')
return_mean_map(res_5_Q,labels[4]+'_residual_Q_mean.fits')
return_mean_map(res_6_Q,labels[5]+'_residual_Q_mean.fits')
return_mean_map(res_7_Q,labels[6]+'_residual_Q_mean.fits')
#return_mean_map(res_8_Q,labels[7]+'_residual_Q_mean.fits')

return_mean_map(res_1_U,labels[0]+'_residual_U_mean.fits')
return_mean_map(res_2_U,labels[1]+'_residual_U_mean.fits')
return_mean_map(res_3_U,labels[2]+'_residual_U_mean.fits')
return_mean_map(res_4_U,labels[3]+'_residual_U_mean.fits')
return_mean_map(res_5_U,labels[4]+'_residual_U_mean.fits')
return_mean_map(res_6_U,labels[5]+'_residual_U_mean.fits')
return_mean_map(res_7_U,labels[6]+'_residual_U_mean.fits')
#return_mean_map(res_8_U,labels[7]+'_residual_U_mean.fits')
