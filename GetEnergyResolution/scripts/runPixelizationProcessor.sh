#!/bin/bash

#not very automatized yet! this is a test only for e- runs from 20240821
source ../../init_ilcsoft.sh

cd /lhome/ific/a/almanzam/testsSim/ExampleProcessor/GetEnergyResolution/scripts
export MARLIN_DLL="$MARLIN_DLL:$PWD/../lib/libGetEnergyResolutionProcessor.so"

declare -i numberEvents=3
AIDAFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/"
AIDAFileName="testECALe_luxe_v0_QGSP_BERT_e-_"

for energy in 1500 
do
    Marlin --global.MaxRecordNumber=${numberEvents} --MyAIDAProcessor.FileName="$AIDAFilePath$AIDAFileName$energy" --global.LCIOInputFiles=/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/PixelizationLCIO/new_ECALe_luxe_v0_QGSP_BERT_e-_1500MeV_0.slcio test.xml
done