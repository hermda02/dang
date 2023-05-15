import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import healpy as hp
import os
from astropy import units as u
from astropy.modeling.models import BlackBody

nside = 16

T_CMB = 2.728*u.K # Fixsen 1996: ApJ 473: 576-587

#T_CMB = 2.72548*u.K # Fixsen 2009, arxiv:0911.1955v2

cmb_bb = BlackBody(temperature=T_CMB)
cmb_map = np.empty(hp.nside2npix(nside))

dirname = 'healpix_dump2'

os.makedirs(dirname,exist_ok=True)
os.makedirs(dirname+'/cmb_monopoles/',exist_ok=True)
os.makedirs(dirname+'/no_monopole/',exist_ok=True)

frequencies = []

fits1 = fits.open('firas_hpx_destriped_sky_spectra_lowf_v2.fits')
fits2 = fits.open('firas_hpx_c_vector_lowf_v2.fits')
nu_zero = float(fits1[0].header['NU_ZERO'])
dnu = float(fits1[0].header['DELTA_NU'])
nfreq = int(fits1[0].header['NUM_FREQ'])
m = np.zeros((3072, nfreq))
r = np.zeros((3072, nfreq))
for i in range(3072):
    for j in range(nfreq):
        m[i,j] = fits1[1].data[i][5][j]
        r[:,j] = fits2[1].data[0][0][j]

for j in range(nfreq):
    freq = float('%.1f'%(nu_zero+j*dnu))
    freq_str =str('%.1f'%(nu_zero+j*dnu)).zfill(6)
    if freq_str not in frequencies:
        frequencies.append(freq_str)
        cmb_map[:] = cmb_bb(freq*u.GHz) 
        hp.write_map(f'{dirname}/cmb_monopoles/cmb_monopole_{freq_str}.fits',cmb_map)
    print(freq_str)
    hp.write_map(f'{dirname}/lowf_map_{freq_str}.fits',m[:,j])
    hp.write_map(f'{dirname}/no_monopole/lowf_map_{freq_str}_no_monopole.fits',m[:,j]-cmb_map)
    hp.write_map(f'{dirname}/lowf_rms_{freq_str}.fits',r[:,j])
    hp.mollview(m[:,j], nest=True, norm='hist')
    plt.savefig(f'{dirname}/lowf_map_{freq_str}.png')
    plt.close()

fits1 = fits.open('firas_hpx_destriped_sky_spectra_high_v2.fits')
fits2 = fits.open('firas_hpx_c_vector_high_v2.fits')
nu_zero = float(fits1[0].header['NU_ZERO'])
dnu = float(fits1[0].header['DELTA_NU'])
nfreq = int(fits1[0].header['NUM_FREQ'])
m = np.zeros((3072, nfreq))
r = np.zeros((3072, nfreq))
for i in range(3072):
    for j in range(nfreq):
        m[i,j] = fits1[1].data[i][5][j]
        r[:,j] = fits2[1].data[0][0][j]

for j in range(nfreq):
    freq = float('%.1f'%(nu_zero+j*dnu))
    freq_str =str('%.1f'%(nu_zero+j*dnu)).zfill(6)
    if freq_str not in frequencies:
        frequencies.append(freq_str)
        cmb_map[:] = cmb_bb(freq*u.GHz) 
        hp.write_map(f'{dirname}/cmb_monopoles/cmb_monopole_{freq_str}.fits',cmb_map)
    print(freq_str)
    hp.write_map(f'{dirname}/highf_map_{freq_str}.fits',m[:,j])
    hp.write_map(f'{dirname}/no_monopole/highf_map_{freq_str}_no_monopole.fits',m[:,j]-cmb_map)
    hp.write_map(f'{dirname}/highf_rms_{freq_str}.fits',r[:,j])
    hp.mollview(m[:,j], nest=True, norm='hist')
    plt.savefig(f'{dirname}/highf_map_{freq_str}.png')
    plt.close()


with open(f'{dirname}/firas_frequencies.txt','w') as f:
    for freq in frequencies:
        f.write("%s\n" % freq)
f.close()
