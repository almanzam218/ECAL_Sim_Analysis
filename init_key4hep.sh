#!/usr/bin/env bash

export GYROS_PATH=$( realpath $(dirname $BASH_SOURCE) )
source $GYROS_PATH/key4hep_latest.env

# TODO: Add the path of LUXE Style
# LUXE_STYLE_PATH=""
# export CPATH="$LUXE_STYLE_PATH:$CPATH"

# TODO: Add all compiled Marlin lib files
# MARLIN_LIB_PATH=""
# for dll in $( ls $MARLIN_LIB_PATH/*.so ); do
#    export MARLIN_DLL="${dll}:$MARLIN_DLL"
# done
