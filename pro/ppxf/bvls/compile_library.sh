#!/bin/bash
# command from https://trac.sdss.org/wiki/MANGA/Data/DAPDevelopment/DAP_v0_9_GettingStarted#Compiletheexternallibraryifdesired
gfortran -dynamiclib -o bvls.dylib  bvls_wrapper.c bvls.f90
