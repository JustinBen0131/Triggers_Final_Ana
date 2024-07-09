In src folder do:

./autogen.sh

./configure

make

make install -j8

to set up using CaloEmulatorTreeMaker.cc

Run Fun4All_Calo_Emulator.C in macros folder:

Go to ‘triggerCondor.sub’

Change run number in arguments and frome queue filename 


Go to ‘triggerCondor.sh’ change MYINSTALL to your install directory 

Make a folder called ‘output’ in directory with ‘macros’, ‘src’ and condor files (remove from condor folder)

Run condor_submit triggerCondor.sub


Go to output/runNumber (folder automatically generated into a path labeled as ‘output’ and see two kinds of files:

CaloOutput_DST_CALO_run2pp_new_2024p001-00044686-*.root, with * = segment number


Go to macros folder where there is createListFile.cpp:
run 'createListFile.cpp' with:
root -b -q -l 'createListFile.cpp("list,of,runs") \\make sure no space between run numbers


This will create a filelist in your base directory with the CaloOutput information to be run in condor.

Can then run ‘CaloTriggerSegments.cpp’ by going to base directory with ‘makeTriggerHistograms.sh’ and ‘makeTriggerHistograms.sub’

You can change any of the filling functionality in ‘CaloTriggerSegments.cpp’ that is needed, the main functionality is it makes histograms of whatever you fill with it, for each trigger index for scaled, and live counts, both with and without scaledown factor applied. 

Edit this functionality as needed.
