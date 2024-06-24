#!/usr/bin/env bash

# Environment setup
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
nEvents=${4:-0}  # Default number of events to 0 if not provided

# Debugging statements
echo "runNumber: $runNumber"
echo "filename: $filename"
echo "clusterID: $clusterID"
echo "nEvents: $nEvents"

# Ensure the filename is a full path
if [[ "$filename" != /* ]]; then
  echo "Error: Filename must be a full path."
  exit 1
fi

# Create a directory based on the run number
outputDir="outputHistRootFiles/${runNumber}"
mkdir -p ${outputDir}  # Ensure the directory exists

# Debugging statements
echo "outputDir: $outputDir"

# Construct the output filenames within the new directory
baseFilename=$(basename -- "$filename")
filenameWithoutExt="${baseFilename%.*}"
treeOutName="${outputDir}/HistogramOutput_${filenameWithoutExt}.root"  # Strip the extension from the filename and append

# Check if filename is provided and exists
if [ -z "$filename" ]; then
  echo "Error: Filename is not provided."
  exit 1
fi

if [ ! -f "$filename" ]; then
  echo "Error: Input file $filename does not exist."
  exit 1
fi

# Debugging statements
echo "Processing rootFile: $filename"
echo "treeOutName: $treeOutName"

# Run the ROOT macro with dynamically generated output paths and log the output
root -b -l -q "macros/CaloTriggerSegments.cpp(\"$filename\", \"$treeOutName\", $nEvents)"
