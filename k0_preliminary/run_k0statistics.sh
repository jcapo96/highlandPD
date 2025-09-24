#!/bin/bash

# Script to run K0Statistics.C analysis
# This script analyzes K0 production from K+ beam particles
# Modified to process all files in DATA directory and produce combined results

echo "Starting K0 Statistics Analysis (Combined from all DATA files)..."
echo "=================================================================="

# Set up ROOT environment
source /Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/highland/scripts/setup.sh

# Manual ROOT setup if automatic setup fails
if ! command -v root >/dev/null 2>&1; then
    echo "ROOT not found in PATH, setting up manually..."
    export ROOTSYS="/Applications/root_src/root_install"
    export PATH="$ROOTSYS/bin:$PATH"
    export DYLD_LIBRARY_PATH="$ROOTSYS/lib/root:$ROOTSYS/lib"
fi

# Change to the correct directory
cd /Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/highlandPD/k0_preliminary

# Check if DATA directory exists
DATA_DIR="/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/DATA"
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: DATA directory not found at $DATA_DIR"
    exit 1
fi

# Count files in DATA directory
file_count=$(ls -1 "$DATA_DIR"/*.root 2>/dev/null | wc -l)
echo "Found $file_count root files in DATA directory"
echo ""

# Create temporary script that loads libraries and calls the function
TEMP_SCRIPT=$(mktemp)
cat > "$TEMP_SCRIPT" << EOF
gSystem->Load("libhighland");
gSystem->Load("libhighlandPD");
.L K0Statistics.C
K0Statistics("$DATA_DIR");
EOF

# Run the analysis
echo "Running K0Statistics.C with combined analysis..."
root -l -b -q "$TEMP_SCRIPT"

# Clean up temporary file
rm -f "$TEMP_SCRIPT"

echo ""
echo "Combined analysis complete!"
echo "Results show statistics from all files processed as a single dataset."
