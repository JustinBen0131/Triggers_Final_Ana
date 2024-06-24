#!/usr/bin/env bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
export MYINSTALL="/sphenix/user/patsfan753/install"

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

# Capture command line arguments
runNumber=$1
filename=$2
clusterID=$3
events=0  # Default number of events

# Create a directory based on the run number
outputDir="output/${runNumber}"
mkdir -p ${outputDir}  # Ensure the directory exists

# Construct the output filenames within the new directory
treeOutName="${outputDir}/CaloOutput_${filename%.*}.root"  # Strip the extension from the filename and append
qaOutputFilename="${outputDir}/QAoutput_${filename%.*}.root"

# Check if filename is provided
if [ -z "$filename" ]; then
  echo "Error: Filename is not provided."
  exit 1
fi

# Run the ROOT macro with dynamically generated output paths and log the output
root -b -l -q "macros/Fun4All_Calo_Emulator.C(\"$filename\", $events, \"$treeOutName\", \"$qaOutputFilename\")"
