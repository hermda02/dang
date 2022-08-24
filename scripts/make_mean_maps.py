import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os
import sys

hp.disable_warnings()

def return_mean_map(list,outfile,mask):
    print('Creating '+str(outfile))
    nside = 64
    npix  = nside*nside*12
    samples = len(list)

    print(samples)

    out_map = np.empty((npix,3))

    out_map[:][:] = 0.0
    
    for i in range(samples):
        out_map = out_map + hp.read_map(list[i],verbose=False,field=(0,1,2)).T

    out_map = out_map/float(samples)
    for i in range(npix):
        if mask[i] == 0:
            for j in range(3):
                out_map[i][j] = hp.UNSEEN 

    hp.write_map(str(outfile),out_map.T)

def return_std_map(list,outfile,mask):
    print('Creating '+str(outfile))
    nside = 64
    npix  = nside*nside*12
    samples = len(list)

    print(samples)

    maps = np.empty((npix,3,samples))
    out_map = np.empty((npix,3))
    
    maps = maps.T
    
    for i in range(samples):
        maps[i] = hp.read_map(list[i],verbose=False,field=(0,1,2))

    maps = maps.T

    for i in range(npix):
        if mask[i] == 0:
            for j in range(3):
                out_map[i][j] = hp.UNSEEN 
        else:
            for j in range(3):
                out_map[i][j] = np.std(maps[i][j][:])


    hp.write_map(outfile,out_map.T)
    
def read_params(filename):
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

try:    
    files = os.listdir()
    thing = [file for file in files if 'param' in file]
    names, freq, num_samp, maskfile = read_params(thing[0])
    labels = [name.replace("'","") for name in names]
    maskfile = maskfile.replace("'","")
    num_bands = len(freq)
    print(labels)
    print(freq)
except:
    print("Input which directory you wish to point to (../ automatically included)")
    exit()

# listed  = [s for s in files if (labels[0] in s and 'res' in s)]
# samples = len(listed)
# print(samples)

mask = hp.read_map(str('../'+maskfile))

# npix = np.shape(mask)[0]

# residuals = np.empty((num_bands,samples,npix))
# synchs    = np.empty((samples,npix))
# chis      = np.empty((samples,npix))
# betas     = np.empty((samples,npix))

# for file in files:
#     for j in range(num_bands):
#         if (labels[j] in file) and ('res' in file):
#             print(labels[j],file)
#     # for i in range(samples):



res_1 = []
res_2 = []
res_3 = []
res_4 = []
res_5 = []
res_6 = []
res_7 = []
# res_8 = []

synchs = []
chis   = []

beta_s  = []

for file in files:
    
    # Residuals Q
    if file.startswith(labels[0]+'_residual'):
        if file.endswith('.fits'):
            res_1.append(file)

    if file.startswith(labels[1]+'_residual'):
        if file.endswith('.fits'):
            res_2.append(file)
    
    if file.startswith(labels[2]+'_residual'):
        if file.endswith('.fits'):
            res_3.append(file)

    if file.startswith(labels[3]+'_residual'):
        if file.endswith('.fits'):
            res_4.append(file)

    if file.startswith(labels[4]+'_residual'):
        if file.endswith('.fits'):
            res_5.append(file)

    if file.startswith(labels[5]+'_residual'):
        if file.endswith('.fits'):
            res_6.append(file)

    if file.startswith(labels[6]+'_residual'):
        if file.endswith('.fits'):
            res_7.append(file)


    # Synch maps
            
    if file.startswith('bp_030_synch_amplitude'):
        if file.endswith('.fits'):
            synchs.append(file)

    # Synch beta
    if file.startswith('synch_beta'):
        if file.endswith('.fits'):
            beta_s.append(file)
            
    # Chisquare maps

    if file.startswith('chisq'):
        if file.endswith('.fits'):
            chis.append(file)


return_mean_map(chis,'chisq_mean.fits',mask)
return_mean_map(beta_s,'synch_beta_mean.fits',mask)
return_mean_map(synchs,'synch_mean.fits',mask)
return_std_map(synchs,'synch_std.fits',mask)

return_std_map(beta_s,'synch_beta_std.fits',mask)

return_mean_map(res_1,labels[0]+'_residual_mean.fits',mask)
return_mean_map(res_2,labels[1]+'_residual_mean.fits',mask)
return_mean_map(res_3,labels[2]+'_residual_mean.fits',mask)
return_mean_map(res_4,labels[3]+'_residual_mean.fits',mask)
return_mean_map(res_5,labels[4]+'_residual_mean.fits',mask)
return_mean_map(res_6,labels[5]+'_residual_mean.fits',mask)
return_mean_map(res_7,labels[6]+'_residual_mean.fits',mask)
# #return_mean_map(res_8,labels[7]+'_residual_mean.fits')
