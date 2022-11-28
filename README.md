# dang

CONTENTS OF THIS FILE
---------------------

* Overview/disclaimer
* Introduction
* Requirements

OVERVIEW/DISCLAIMER
-------------------
Gibbs sampler created for microwave sky component separation.

This version has undergone a massive restructuring to prepare the code for future usability and smarter object orientation.

The current head of the `master` branch reproduces the results of arxiv:2201.03530, but the original paper was published using the following hash:
512fbd7fdbde93ef45b704a8edbebc57fd72256e

INTRODUCTION
------------
This code has been developed for Bayesian data analysis of the microwave/sub-mm sky in pixel space. 
This code utilizes Gibbs sampling to map out the distribution of parameters. Functionalities are 
regularly added as needed.


Requirements
------------

This code has been developed in FORTRAN 90, despite the ancient nature of the language.

Requirements are as follows:

- CFITSIO (https://heasarc.gsfc.nasa.gov/fitsio/)
- HEALPix (and dependencies) (https://healpix.sourceforge.io/) 
- LAPACK/BLAS
- OpenMP

If one wishes to utilize the full speed of this software, they must manually link parallelized LAPACK
libraries within the Makefile.