#!/bin/bash

# Script to run Vertex Display Analysis
# Combines event iteration with vertex creation algorithm to visualize events

# Set up the environment
echo "Setting up HIGHLAND environment..."

# Source the setup script - use absolute path like working scripts
echo "Sourcing HIGHLAND setup.sh..."
source /Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/highland/scripts/setup.sh

# Manual ROOT setup if automatic setup fails (like in working scripts)
if ! command -v root >/dev/null 2>&1; then
    echo "ROOT not found in PATH, setting up manually..."
    export ROOTSYS="/Applications/root_src/root_install"
    export PATH="$ROOTSYS/bin:$PATH"
    export DYLD_LIBRARY_PATH="$ROOTSYS/lib/root:$ROOTSYS/lib"
fi

# Change to the correct directory (like working scripts)
cd /Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/highlandPD/k0_preliminary

echo ""
echo "Environment setup complete."
echo "ROOTSYS: $ROOTSYS"
echo "HIGHLANDROOT: $HIGHLANDROOT"
echo ""

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <input_file> [max_events]"
    echo "Example: $0 /path/to/minitree.root 5"
    echo "Default max_events: 10"
    exit 1
fi

INPUT_FILE="$1"
MAX_EVENTS="${2:-1000}"

echo "Input file: $INPUT_FILE"
echo "Maximum events to process: $MAX_EVENTS"
echo ""

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

# Run the analysis
echo "Running Vertex Display Analysis..."
echo "This will create visualization images for events with vertices."
echo ""

# Create temporary script that loads libraries and calls the function (like working scripts)
TEMP_SCRIPT=$(mktemp)
cat > "$TEMP_SCRIPT" << EOF
gSystem->Load("libhighland");
gSystem->Load("libhighlandPD");
.L VertexDisplayAnalysis.C
VertexDisplayAnalysis("$INPUT_FILE", $MAX_EVENTS);
EOF

# Run ROOT with the script
root -l -b -q "$TEMP_SCRIPT"

# Clean up temporary file
rm -f "$TEMP_SCRIPT"

echo ""
echo "Analysis complete!"
echo "Check for generated visualization files: event_*_visualization.png"
