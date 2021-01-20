# spfit

## Introduction

SPFIT is an IDL package to simultaneous fit stellar continuum and emission lines. The package is built upon the Penalized Pixel-Fitting method (pPXF) (Cappellari 2017) and Gas AND Absorption Line Fitting (GANDALF) (Sarzi et al. 2006). For the ease of use, wrappers and examples for modeling SDSS single-fiber spectra and MaNGA IFU datacubes are provided. 

## References

Fu et al. 2018, https://ui.adsabs.harvard.edu/abs/2018ApJ...856...93F

Steffen et al. 2021, ApJ, currently under review. 

## Setup

1. Download SPFIT 
```cs
cd ~/idl
git clone https://github.com/fuhaiastro/spfit.git
```
2. Setup C shell environment
add the following commend. This is used to setup correct environment variables when launching IDL *
```cs
vi ~/.tcshrc
alias spfit 'source ~/idl/spfit/spfit.csh; /Applications/harris/idl/bin/idl -IDL_PROMPT "SPFIT> " -IDL_STARTUP ""'
```
3. Edit ~/idl/spfit/spfit.csh 

Make sure the appropriate paths are set correctly for your host.
```cs
vi ~/idl/spfit/spfit.csh
```
4. INSTALL IDLUTILS 
```cs
cd ~/idl
svn co https://www.sdss3.org/svn/repo/idlutils/trunk/ idlutils
setenv IDLUTILS_DIR $IDL/idlutils
$IDLUTILS_DIR/bin/evilmake clean
$IDLUTILS_DIR/bin/evilmake 
```
5. INSTALL IDLSPEC2D
```cs
cd ~/idl
svn co https://www.sdss3.org/svn/repo/idlspec2d/trunk/ idlspec2d
setenv IDLSPEC2D_DIR $IDL/idlspec2d
$IDLSPEC2D_DIR/bin/evilmake clean
$IDLSPEC2D_DIR/bin/evilmake 
```

## Examples

Examples for SDSS/BOSS and MaNGA are provided under `spfit/examples/`
See the README files in these folders for instructions.

