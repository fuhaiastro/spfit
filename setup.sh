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
