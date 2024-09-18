#!/bin/bash

#not very automatized yet! this is a test only for e- runs from 20240821
source ../../init_ilcsoft.sh

cd /lhome/ific/a/almanzam/testsSim/ExampleProcessor/GetPIDParametersTree/scripts
export MARLIN_DLL="$MARLIN_DLL:$PWD/../lib/libGetPIDParametersTreeProcessor.so"

declare -i numberEvents=3

for energy in 10 
do
    AIDAFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/PIDParametersTrees/tests/20240917/"
    AIDAFileName="PIDParams_PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_gamma_${energy}GeV"

    InputLCIOFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/NPOD_samples/20240829_v1/clustering/data/"
    InputLCIOFileName="PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_gamma_${energy}GeV.slcio"
    Marlin --global.MaxRecordNumber=${numberEvents} --MyAIDAProcessor.FileName="$AIDAFilePath$AIDAFileName" --global.LCIOInputFiles=$InputLCIOFilePath$InputLCIOFileName test.xml
done