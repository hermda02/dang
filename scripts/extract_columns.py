import healpy as hp
import numpy as np
import sys
import os

from astropy.io import fits

def dump_TQU(data,nfields,unit):
    if nfields == 10:

        if unit == 'Kcmb':
            I_map = data[1].data['I_STOKES']
            Q_map = data[1].data['Q_STOKES']
            U_map = data[1].data['U_STOKES']
            
            I_cov = data[1].data['II_COV']
            Q_cov = data[1].data['QQ_COV']
            U_cov = data[1].data['UU_COV']

        # Missing pixel check
        npix = len(I_map)
        
        map_inds = (I_map == hp.UNSEEN)
        cov_inds = (I_cov == 0.0)

        I_map *= 1e6
        Q_map *= 1e6
        U_map *= 1e6

        I_map[map_inds] = hp.UNSEEN
        Q_map[map_inds] = hp.UNSEEN
        U_map[map_inds] = hp.UNSEEN

        # hp.reorder(I_map,inp='NESTED',out='RING')
        # hp.reorder(Q_map,inp='NESTED',out='RING')
        # hp.reorder(U_map,inp='NESTED',out='RING')

        # hp.reorder(I_cov,n2r=True)
        # hp.reorder(Q_cov,n2r=True)
        # hp.reorder(U_cov,n2r=True)

        hp.write_map("HFI_map_"+freq+"_IQU_n2048_u"+unit+".fits",(I_map,Q_map,U_map),overwrite=True,nest=True)
        hp.write_map("HFI_rms_"+freq+"_IQU_n2048_u"+unit+".fits",(np.sqrt(I_cov)*1e6,np.sqrt(Q_cov)*1e6,np.sqrt(U_cov)*1e6),overwrite=True,nest=True)

    elif nfields == 3:
        unit = str(unit.split('/')[0])
        
        I_map = data[1].data['I_STOKES']
        I_cov = data[1].data['II_COV']
        
        hp.write_map("HFI_map_"+freq+"_I_n2048_"+unit+".fits",(I_map),overwrite=True,nest=True)
        hp.write_map("HFI_rms_"+freq+"_I_n2048_"+unit+".fits",(np.sqrt(I_cov)),overwrite=True,nest=True)

    else:
        print(f"That ain't gonna work, nfields = {nfields}.")


files = os.listdir()

skymaps = [file for file in files if 'SkyMap' in file]
skymaps.sort()

for mapname in skymaps:

    print(mapname)
    
    freq = str(mapname.split('_')[2])
    hdul = fits.open(mapname)
    hdr = hdul[1].header
    
    n_fields  = hdr['TFIELDS']
    unit = hdr['TUNIT1']

    dump_TQU(hdul,n_fields,unit)

    print('')
    # exit()
