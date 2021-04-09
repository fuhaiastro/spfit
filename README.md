# spfit

## Introduction

SPFIT is an IDL package to simultaneous fit stellar continuum and emission lines. The package is built upon the Penalized Pixel-Fitting method (pPXF) (Cappellari 2017) and Gas AND Absorption Line Fitting (GANDALF) (Sarzi et al. 2006). For the ease of use, wrappers and examples for modeling SDSS single-fiber spectra and MaNGA IFU datacubes are provided.

## References

Fu et al. 2018, https://ui.adsabs.harvard.edu/abs/2018ApJ...856...93F \
Steffen et al. 2021, https://arxiv.org/abs/2102.03398

## Setup

1. Download SPFIT
```shell
cd ~/idl
git clone https://github.com/fuhaiastro/spfit.git
```
2. Setup C shell environment
add the following commend. This is used to setup correct environment variables when launching IDL *
```shell
vi ~/.tcshrc
alias spfit 'source ~/idl/spfit/setup.csh; /Applications/harris/idl/bin/idl -IDL_PROMPT "SPFIT> "'
```
3. Edit ~/idl/spfit/spfit.csh

Make sure the appropriate paths are set correctly for your host.
```shell
vi ~/idl/spfit/spfit.csh

# IDL environment variables and aliases.
# tested on IDL v8.6 for Mac OSX
setenv EXELIS_DIR /Applications/harris
setenv IDL_DIR /Applications/harris/idl
alias harrislicense $IDL_DIR/bin/harrislicense
if ( -x $IDL_DIR/bin/idlde ) alias idlde $IDL_DIR/bin/idlde
if ( -x $IDL_DIR/bin/idlhelp ) alias idlhelp $IDL_DIR/bin/idlhelp
if ( -x $IDL_DIR/bin/idlrpc ) alias idlrpc $IDL_DIR/bin/idlrpc
if ( -x $IDL_DIR/bin/idltaskengine ) alias idltaskengine $IDL_DIR/bin/idltaskengine
setenv IDL_PATH +$IDL_DIR/examples:+$IDL_DIR/lib
setenv IDL $HOME/idl

# MaNGA
setenv MANGA_DIR /s1/manga # Data Directory
setenv MANGA_SPECTRO_REDUX $MANGA_DIR/spectro/redux/
setenv MANGADRP_VER MPL-11
# SPFIT
setenv SPFIT_DIR $HOME/idl/spfit/
setenv IDL_PATH ${IDL_PATH}:+$SPFIT_DIR/pro
# IDLUTILS (includes Astrolib,Coyote,MPFIT)
setenv IDLUTILS_DIR $IDL/idlutils
setenv DUST_DIR $IDL/dust # SFD98 maps
setenv IDL_PATH ${IDL_PATH}:+$IDLUTILS_DIR/goddard/pro:+$IDLUTILS_DIR/pro
setenv PATH $IDLUTILS_DIR/bin:$PATH
# IDLSPEC2D
setenv IDLSPEC2D_DIR $IDL/idlspec2d
setenv IDL_PATH ${IDL_PATH}:+$IDLSPEC2D_DIR/pro
setenv PATH $IDLSPEC2D_DIR/bin:$PATH
```
4. Install IDLUTILS
```shell
cd ~/idl
mkdir idlutils
cd idlutils
svn co https://svn.sdss.org/public/repo/sdss/idlutils/trunk .
setenv IDL_DIR /Applications/harris/idl
setenv IDLUTILS_DIR $HOME/idl/idlutils
$IDLUTILS_DIR/bin/evilmake clean
$IDLUTILS_DIR/bin/evilmake
```
5. Install IDLSPEC2D
```shell
cd ~/idl
mkdir idlspec2d
cd idlspec2d
svn co https://svn.sdss.org/public/repo/eboss/idlspec2d/trunk/ .
setenv IDLSPEC2D_DIR $IDL/idlspec2d
$IDLUTILS_DIR/bin/evilmake clean
$IDLUTILS_DIR/bin/evilmake
```

## Examples

Examples are provided for fitting [SDSS single-fiber
spectra](https://github.com/fuhaiastro/spfit/tree/main/examples/sdss)
and [MaNGA
datacubes](https://github.com/fuhaiastro/spfit/tree/main/examples/manga).
See the README files in these folders for instructions.
