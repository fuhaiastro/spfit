# spfit: a spectral fitting tool in IDL

SPFIT is an IDL-based package to simultaneous fit stellar continuum and emission lines. The package is built upon the Penalized Pixel-Fitting method (pPXF) (Cappellari & Emsellem 2004, Cappellari 2017) and Gas AND Absorption Line Fitting (GANDALF) (Sarzi et al. 2006). For the ease of use, wrappers and examples for modeling SDSS single-fiber spectra and MaNGA IFU datacubes are provided.

## References

- Cappellari & Emsellem [2004PASP..116..138C](http://adsabs.harvard.edu/abs/2004PASP..116..138C)
- Sarzi et al. [2006MNRAS.366.1151S](https://ui.adsabs.harvard.edu/abs/2006MNRAS.366.1151S)
- Cappellari [2017MNRAS.466..798C](https://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C)
- Fu et al. [2018ApJ...856...93F](https://ui.adsabs.harvard.edu/abs/2018ApJ...856...93F)
- Steffen et al. [2021ApJ...909..120S](https://ui.adsabs.harvard.edu/abs/2021ApJ...909..120S)

## Install

1. Download source code
```shell
cd ~/idl
git clone https://github.com/fuhaiastro/spfit.git
```
2. Set environmental variables and define start-up command to launch `spfit`.
```shell
vi ~/.tcshrc
alias spfit 'source ~/idl/spfit/setup.csh; /Applications/harris/idl/bin/idl -IDL_PROMPT "SPFIT> "'
```
3. Make sure the appropriate paths are set correctly for your host.
```shell
vi ~/idl/spfit/setup.csh
```
Below is the content of the setup file:
```shell
# IDL environment variables and aliases.
# tested on IDL v8.7.2 for Mac OSX
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
setenv MANGA_DIR /s1/manga
setenv MANGA_SPECTRO_REDUX $MANGA_DIR/spectro/redux/
setenv MANGADRP_VER MPL-11
# SPFIT
setenv SPFIT_DIR $HOME/idl/spfit/
setenv IDL_PATH ${IDL_PATH}:+$SPFIT_DIR/pro
# IDLUTILS (includes Astrolib, Coyote, & MPFIT)
setenv IDLUTILS_DIR $IDL/idlutils
setenv IDL_PATH ${IDL_PATH}:+$IDLUTILS_DIR/goddard/pro:+$IDLUTILS_DIR/pro
setenv PATH $IDLUTILS_DIR/bin:$PATH
# SFD98 Galactic dust maps
setenv DUST_DIR $IDL/dust
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
setenv IDL $HOME/idl
setenv IDLUTILS_DIR $IDL/idlutils
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
