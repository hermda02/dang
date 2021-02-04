import healpy as hp
import numpy as np 
import os
import sys

try:    
    dir = sys.argv[1] 
    print(dir)
except:
    print("Input which directory you wish to point to (../ automatically included)")
    exit()

def write_diff(title,map):
    hp.write_map('../'+dir+'/'+title,map)

scale_to_30 = (30./44.)**(-3.1)
scale_to_spass = (2.305/44.)**(-3.1)

joint_synch_Q = hp.read_map('../'+dir+'/synch_Q_mean.fits')
joint_synch_U = hp.read_map('../'+dir+'/synch_U_mean.fits')

joint_Q_30 = scale_to_30*joint_synch_Q
joint_U_30 = scale_to_30*joint_synch_U
joint_Q_spass = scale_to_spass*joint_synch_Q
joint_U_spass = scale_to_spass*joint_synch_U

bp_030_Q = hp.read_map('../data/BP_synch_Q_n0064.fits')
bp_030_U = hp.read_map('../data/BP_synch_U_n0064.fits')

npipe_30_Q = hp.read_map('../data/npipe6v20_comm_synch_n0064_60arc_Q_rc1.fits')
npipe_30_U = hp.read_map('../data/npipe6v20_comm_synch_n0064_60arc_U_rc1.fits')

spass_Q    = hp.read_map('../data/spass_rmrot_n0064_ring_masked.fits',field=1)
spass_U    = hp.read_map('../data/spass_rmrot_n0064_ring_masked.fits',field=2)

bp_min_joint_Q = bp_030_Q - joint_Q_30
bp_min_joint_U = bp_030_U - joint_U_30

write_diff('BP_minus_joint_Q_60arcmin_n0064.fits',bp_min_joint_Q)
write_diff('BP_minus_joint_U_60arcmin_n0064.fits',bp_min_joint_U)

np_min_joint_Q = npipe_30_Q - joint_Q_30
np_min_joint_U = npipe_30_U - joint_U_30

write_diff('npipe_minus_joint_Q_60arcmin_n0064.fits',np_min_joint_Q)
write_diff('npipe_minus_joint_U_60arcmin_n0064.fits',np_min_joint_U)

spass_min_joint_Q = spass_Q - joint_Q_spass
spass_min_joint_U = spass_U - joint_U_spass

write_diff('spass_minus_joint_Q_60arcmin_n0064.fits',spass_min_joint_Q)
write_diff('spass_minus_joint_U_60arcmin_n0064.fits',spass_min_joint_U)
