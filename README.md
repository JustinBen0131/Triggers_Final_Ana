Copy folder from Dans repository:
https://github.com/sPHENIX-Collaboration/analysis/blob/a8c0d512e275232293a9008635aedf2826fc193f/calotriggeremulator


cd into src folder and do:
sed -i '37s/useLL1=false;/useLL1=true;/' /path/to/folder/calotriggeremulator/src/CaloEmulatorTreeMaker.cc

which changes: useLL1=true; → line 37

Run: 
./autogen.sh

./configure

make

make install -j8



Copy this repositories condor files, and macros folder (replacing with current 'macros')

To run Fun4All_Calo_Emulator.C in macros folder:

Go to ‘triggerCondor.sub’

Change run number in arguments and frome queue filename 


Go to ‘triggerCondor.sh’ change MYINSTALL to your install directory 

Make a folder called ‘output’ in directory with ‘macros’, ‘src’ and condor files

Go into triggerCondor_submit.sh and change the list of run numbers for what ever dst lists you have in your directory


Output generates to output/runNumber:

CaloOutput_DST_CALO_run2pp_new_2024p001-00044686-*.root, with * = segment number (equal to #dst root files in calo dst list)


Go to macros folder where there is createListFile.cpp:
run 'createListFile.cpp' with:
root -b -q -l 'createListFile.cpp("list,of,runs") \\make sure no space between run numbers

This will create a filelist in your base directory with the CaloOutput information to be run in condor.

Can then run ‘CaloTriggerSegments.cpp’ by going to base directory with ‘makeTriggerHistograms.sh’ via:
go into makeTriggerHistograms_submitJobs.sh and change list of run numbers for Calo Dsts generated in previous jobs

You can change any of the filling functionality in ‘CaloTriggerSegments.cpp’ that is needed 
