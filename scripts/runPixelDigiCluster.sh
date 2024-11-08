#!/bin/bash

#not very automatized yet! this is a test only for e- runs from 20240821
source ../../init_ilcsoft.sh

cd /lhome/ific/a/almanzam/testsSim/ECALe-lcio/PixelizationProcessor/scripts
export MARLIN_DLL="$MARLIN_DLL:$PWD/../lib/libPixelizationProcessor.so"

declare -i numberEvents=10001

for energy in 1500 3500 5500 7500 9500 11500 13500 15000
do
    PixelatedFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/PixelizationLCIO/"
    PixelatedFileName="new_ECALe_luxe_v0_QGSP_BERT_e-_${energy}MeV_0.slcio"

    InputLCIOFilePath="/lhome/ific/a/almanzam/gluon/LUXE/ddsim/electron/run_20240821/"
    InputLCIOFileName="ECALe_luxe_v0_QGSP_BERT_e-_${energy}MeV_0.slcio"
    Marlin --global.MaxRecordNumber=${numberEvents} --DSTOutput.LCIOOutputFile="$PixelatedFilePath$PixelatedFileName" --global.LCIOInputFiles=$InputLCIOFilePath$InputLCIOFileName test.xml
done
