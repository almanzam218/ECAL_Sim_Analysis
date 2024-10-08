#!/bin/bash

#not very automatized yet! this is a test only for e- runs from 20240821
source ../../init_ilcsoft.sh

cd /lhome/ific/a/almanzam/testsSim/ExampleProcessor/GetPIDParametersTree/scripts
export MARLIN_DLL="$MARLIN_DLL:$PWD/../lib/libGetPIDParametersTreeProcessor.so"

declare -i numberEvents=2
 #--global.MaxRecordNumber=${numberEvents}
for energy in 10 
do
    particleName="neutron"
    AIDAFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/PIDParametersTrees/20240827_v1/"
    AIDAFileName="PIDParams_PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_${particleName}_0.5to10GeV"

    InputLCIOFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/NPOD_samples/20240829_v1/clustering/data/"
    InputLCIOFileName="PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_${particleName}_${energy}GeV.slcio"
    Marlin  --MyAIDAProcessor.FileName="$AIDAFilePath$AIDAFileName"  test.xml
done