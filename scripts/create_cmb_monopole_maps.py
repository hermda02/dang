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

cmb_spectrum = cmb_bb(frequencies*u.GHz).to('MJy/sr')

cmb_monopole = np.empty(hp.nside2npix(nside))
for i in range(len(cmb_spectrum)):
    cmb_monopole[:] = cmb_spectrum[i]
    hp.write_map(f'cmb_monopole_{frequencies[i]}.fits',cmb_monopole)

