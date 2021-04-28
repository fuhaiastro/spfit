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
