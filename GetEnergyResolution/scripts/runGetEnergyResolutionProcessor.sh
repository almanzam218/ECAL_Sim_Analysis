#!/bin/bash

#not very automatized yet! this is a test only for e- runs from 20240821
source ../../init_ilcsoft.sh

cd /lhome/ific/a/almanzam/testsSim/ExampleProcessor/GetEnergyResolution/scripts
export MARLIN_DLL="$MARLIN_DLL:$PWD/../lib/libGetEnergyResolutionProcessor.so"

declare -i numberEvents=10001

for energy in 1500 3500 5500 7500 9500 11500 13500 15000 #1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000 10500 11000 11500 12000 12500 13000 13500 14000 14500 15000
do
    AIDAFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/20240827_v1_pixelised/"
    AIDAFileName="Pixelization_ECALe_luxe_v1_QGSP_BERT_e-_${energy}MeV_0"

    InputLCIOFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/NPOD_samples/20240827_v1/EnergyR/pixelization/data/"
    InputLCIOFileName="Pixelization_ECALe_luxe_v1_QGSP_BERT_e-_${energy}MeV_0.slcio"
    Marlin --global.MaxRecordNumber=${numberEvents} --MyAIDAProcessor.FileName="$AIDAFilePath$AIDAFileName" --global.LCIOInputFiles=$InputLCIOFilePath$InputLCIOFileName test.xml
done