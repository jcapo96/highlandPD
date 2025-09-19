#!/bin/bash

# Script to run K0Statistics.C analysis
# This script analyzes K0 production from K+ beam particles

echo "Starting K0 Statistics Analysis..."
echo "=================================="

# Set up ROOT environment
source /Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/highland/scripts/setup.sh

# Change to the correct directory
cd /Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/highlandPD/k0_preliminary

# Run the analysis
echo "Running K0Statistics.C..."
root -l -b -q K0Statistics.C

echo "Analysis complete!"
