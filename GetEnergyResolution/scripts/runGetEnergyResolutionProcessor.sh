#!/bin/bash

#not very automatized yet! this is a test only for e- runs from 20240821
source ../../init_ilcsoft.sh

cd /lhome/ific/a/almanzam/testsSim/ExampleProcessor/GetEnergyResolution/scripts
export MARLIN_DLL="$MARLIN_DLL:$PWD/../lib/libGetEnergyResolutionProcessor.so"

declare -i numberEvents=4

for energy in 1500 2000 
do
    AIDAFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/"
    AIDAFileName="test1ECALe_luxe_v0_QGSP_BERT_e-_${energy}"

    InputLCIOFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/PixelizationLCIO/"
    InputLCIOFileName="new_ECALe_luxe_v0_QGSP_BERT_e-_${energy}MeV_0.slcio"
    Marlin --global.MaxRecordNumber=${numberEvents} --MyAIDAProcessor.FileName="$AIDAFilePath$AIDAFileName" --global.LCIOInputFiles=$InputLCIOFilePath$InputLCIOFileName test.xml
done