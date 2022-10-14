import healpy as hp
import numpy as np

chisq_q = hp.read_map('chisq_mean.fits',field=1)
chisq_u = hp.read_map('chisq_mean.fits',field=2)

chisq_av = (chisq_q+chisq_u)/2.0

hp.write_map('chisq_mean_average.fits',chisq_av)
