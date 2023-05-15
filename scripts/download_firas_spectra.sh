#!/bin/bash

wget=/usr/bin/wget

lowspec_map=https://lambda.gsfc.nasa.gov/data/cobe/hpx_firas/destriped_spectra/firas_hpx_destriped_sky_spectra_lowf_v2.fits
highspec_map=https://lambda.gsfc.nasa.gov/data/cobe/hpx_firas/destriped_spectra/firas_hpx_destriped_sky_spectra_high_v2.fits

lowspec_rms=https://lambda.gsfc.nasa.gov/data/cobe/hpx_firas/c_vector/firas_hpx_c_vector_lowf_v2.fits
highspec_rms=https://lambda.gsfc.nasa.gov/data/cobe/hpx_firas/c_vector/firas_hpx_c_vector_high_v2.fits

$wget $lowspec_map
$wget $highspec_map

$wget $lowspec_rms
$wget $highspec_rms
