#!/bin/bash

#not very automatized yet! this is a test only for e- runs from 20240821
source ../../init_ilcsoft.sh

cd /lhome/ific/a/almanzam/testsSim/ExampleProcessor/GetPIDParametersTree/scripts
export MARLIN_DLL="$MARLIN_DLL:$PWD/../lib/libGetPIDParametersTreeProcessor.so"

declare -i numberEvents=2
 #--global.MaxRecordNumber=${numberEvents}
for energy in 10 
do
    particleName="pi-"
    AIDAFilePath="/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/PIDParametersTrees/20241106_v1/EventDisplay_pi-_20241112/"
    AIDAFileName="PIDParams_PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_${particleName}_Lin_0.5-3GeV_50-54"

    Marlin  --MyAIDAProcessor.FileName="$AIDAFilePath$AIDAFileName"  test.xml
done