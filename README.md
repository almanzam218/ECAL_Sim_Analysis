# Example Processor
Simple processor that reads LCIO files with SimCalorimeterHit and writes in screen the energy deposition and coordinates of each hit

# Compile:

> cd XXXXX/ExampleProcessor
>
> source init_ilcsoft.sh

(or use an updated one from /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03 or similar folders!)

> mkdir build
> 
> cd build
> 
> cmake -C $ILCSOFT/ILCSoft.cmake ..
> 
> make -j3
> 
> make install
>changes to readme

# Run: 
scripts_condor are copied from other example, not usable for us unless few modifications are done.

> cd scripts/
> 
> export MARLIN_DLL="$MARLIN_DLL:$PWD/../lib/libExampleProcessor.so"
> 
(the export MARLIN_DLL is done only once per session, same as the source init_ilcsoft.sh)
> Marlin test.xml
> 
# ECAL_Sim_Analysis
