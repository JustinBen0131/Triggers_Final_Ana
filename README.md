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

And 

QAoutput_DST_CALO_run2pp_new_2024p001-00044686-*.root where * = segment number

First move all the ‘QAoutput..’ into a different folder:

cd output/runnumber

mkdir qaOutput

mv QAoutput_DST_CALO_run2pp_new_2024p001-00044686-*.root .


And then go to macros folder where there is createListFile.cpp

Compile the program:

g++ -std=c++17 -o createListFile createListFile.cpp



Enter:

./createListFile /sphenix/u/patsfan753/scratch/analysis/calotriggeremulator/output/44630 44630



This will create a filelist in your base directory with the CaloOutput information to be run in condor.

Can then run ‘CaloTriggerSegments.cpp’ by going to base directory with ‘makeTriggerHistograms.sh’ and ‘makeTriggerHistograms.sub’

You can change any of the filling functionality in ‘CaloTriggerSegments.cpp’ that is needed, the main functionality is it makes histograms of whatever you fill with it, for each trigger index for scaled, and live counts, both with and without scaledown factor applied. 

Edit this functionality as needed.
