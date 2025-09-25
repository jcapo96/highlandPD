#!/bin/bash

# Script to run start/end inversion analysis
# Analyzes how often reconstructed particles have inverted start/end positions compared to true particles

# Set the data directory (can be overridden by command line argument)
DATA_DIR="/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/DATA"

# Check if a different data directory was provided
if [ $# -gt 0 ]; then
    DATA_DIR="$1"
fi

echo "Running Start/End Inversion Analysis..."
echo "Data directory: $DATA_DIR"
echo ""

# Check if the data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory '$DATA_DIR' does not exist!"
    exit 1
fi

# Check if there are any .root files in the directory
if [ -z "$(find "$DATA_DIR" -name "*.root" -type f)" ]; then
    echo "Error: No .root files found in '$DATA_DIR'!"
    exit 1
fi

# Set up ROOT environment
echo "Setting up ROOT environment..."

# Source the setup script if it exists
if [ -f "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/setup.sh" ]; then
    source /Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/setup.sh
fi

# Run the analysis
echo "Starting analysis..."
echo "=========================================="

root -l -b -q "StartEndInversionAnalysis.C(\"$DATA_DIR\")"

echo ""
echo "Analysis completed!"
echo "=========================================="
