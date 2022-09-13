# spfit: a spectral fitting tool for galaxies

SPFIT is an
[IDL](https://www.l3harrisgeospatial.com/Software-Technology/IDL)-based package to simultaneous fit models of stellar
continuum and emission lines to galaxy spectra. The
intial solution is found with the Penalized Pixel-Fitting method
([pPXF](https://www-astro.physics.ox.ac.uk/~cappellari/software/#ppxf)) and the optimization of final solution uses [MPFIT - Robust non-linear least squares curve
fitting](https://pages.physics.wisc.edu/~craigm/idl/fitting.html). For
the ease of use, wrappers and examples for modeling SDSS single-fiber
spectra and MaNGA integral-field spectroscopic datacubes are provided. 
A brief installation guide can be found below, which assumes macOS and bash shell.

## References

- Steffen et al. [2021ApJ...909..120S](https://ui.adsabs.harvard.edu/abs/2021ApJ...909..120S)
- Fu et al. [2018ApJ...856...93F](https://ui.adsabs.harvard.edu/abs/2018ApJ...856...93F)

## Install

1. Download source code
```shell
cd ~/idl
git clone https://github.com/fuhaiastro/spfit.git
```
2. Set environmental variables and define start-up command to launch `spfit`.
```shell
vi ~/.bash_profile
alias spfit='source ~/idl/spfit/setup.sh; /Applications/harris/idl/bin/idl -IDL_PROMPT "SPFIT>"'
```
3. Make sure the appropriate paths are set correctly for your host.
```shell
vi ~/idl/spfit/setup.sh
```
Below is the content of the setup file:
```shell
# Bash shell commands to define IDL environment variables and aliases.
. /Applications/harris/idl/bin/idl_setup.bash

export IDL_PATH=+$IDL_DIR/examples:+$IDL_DIR/lib
export IDL=$HOME/idl

# MaNGA 
export MANGA_DIR=/Volumes/scr/manga
export MANGA_SPECTRO_REDUX=$MANGA_DIR/spectro/redux/
export MANGADRP_VER=MPL-11
# SPFIT
export SPFIT_DIR=$HOME/idl/spfit/
export IDL_PATH=${IDL_PATH}:+$SPFIT_DIR/pro
# IDLUTILS (includes Astrolib, Coyote, & MPFIT)
export IDLUTILS_DIR=$IDL/idlutils
#export IDL_PATH ${IDL_PATH}:+$IDL/astron/pro:+$IDLUTILS_DIR/goddard/pro:+$IDLUTILS_DIR/pro
export IDL_PATH=${IDL_PATH}:+$IDLUTILS_DIR/goddard/pro:+$IDLUTILS_DIR/pro
export PATH=$IDLUTILS_DIR/bin:$PATH
# SFD98 Galactic dust maps
export DUST_DIR=$IDL/dust 
# IDLSPEC2D
export IDLSPEC2D_DIR=$IDL/idlspec2d
export IDL_PATH=${IDL_PATH}:+$IDLSPEC2D_DIR/pro
export PATH=$IDLSPEC2D_DIR/bin:$PATH
```
4. Install IDLUTILS
```shell
cd ~/idl
mkdir idlutils
cd idlutils
svn co https://svn.sdss.org/public/repo/sdss/idlutils/trunk .
export IDL_DIR=/Applications/harris/idl
export IDL=$HOME/idl
export IDLUTILS_DIR=$IDL/idlutils
$IDLUTILS_DIR/bin/evilmake clean
$IDLUTILS_DIR/bin/evilmake
```
5. Install IDLSPEC2D
```shell
cd ~/idl
mkdir idlspec2d
cd idlspec2d
svn co https://svn.sdss.org/public/repo/eboss/idlspec2d/trunk/ .
export IDLSPEC2D_DIR=$IDL/idlspec2d
$IDLUTILS_DIR/bin/evilmake clean
$IDLUTILS_DIR/bin/evilmake
```
6. Install DUST (Galactic dust maps)
```shell
cd ~/idl
mkdir dust
cd dust
svn co https://svn.sdss.org/public/data/sdss/catalogs/dust/trunk/ .
```

## Demo Examples

Examples are provided for fitting [SDSS single-fiber
spectra](https://github.com/fuhaiastro/spfit/tree/main/examples/sdss)
and [MaNGA
datacubes](https://github.com/fuhaiastro/spfit/tree/main/examples/manga).
See the README files in these folders for instructions.
