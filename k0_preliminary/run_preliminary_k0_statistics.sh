#!/bin/bash

# Script to run Preliminary K0 Selection Statistics
# This script analyzes the efficiency and purity of the K0 selection criteria

# Default file
DEFAULT_FILE="/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root"

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_FILE="$SCRIPT_DIR/PreliminaryK0SelectionStatistic.C"

# Check if the script file exists
if [[ ! -f "$SCRIPT_FILE" ]]; then
    echo "Error: $SCRIPT_FILE not found!"
    exit 1
fi

# Handle command line arguments
if [[ $# -eq 0 ]]; then
    FILE="$DEFAULT_FILE"
elif [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 [filename]"
    echo "  filename    Input ROOT file (default: $DEFAULT_FILE)"
    echo "  -h, --help  Show this help message"
    echo ""
    echo "This script performs statistical analysis of the K0 selection criteria:"
    echo "- Two daughters of the beam particle whose origin is 10-50cm away from the end of the beam particle"
    echo "- These two daughters start positions are closer than 1cm"
    echo "- Checks for K0 in truth information that is a daughter of the beam particle"
    echo ""
    echo "The analysis provides:"
    echo "- Total event counts at each selection stage"
    echo "- K0 efficiency (events with K0 that pass criteria)"
    echo "- Purity (events passing criteria that have K0)"
    echo "- Lists of events with K0 in truth and events passing criteria"
    exit 0
else
    FILE="$1"
fi

# Check if the input file exists
if [[ ! -f "$FILE" ]]; then
    echo "Error: Input file $FILE not found!"
    exit 1
fi

echo "Running Preliminary K0 Selection Statistics..."
echo "Input file: $FILE"
echo "Script: $SCRIPT_FILE"
echo ""

# Create a temporary ROOT macro to load libraries and run the script
TEMP_MACRO=$(mktemp)
cat > "$TEMP_MACRO" << 'EOF'
gSystem->Load("libhighland");
gSystem->Load("libhighlandPD");
.L PreliminaryK0SelectionStatistic.C
PreliminaryK0SelectionStatistic("INPUT_FILE");
.q
EOF

# Replace INPUT_FILE placeholder with actual file path
sed -i.bak "s|INPUT_FILE|$FILE|g" "$TEMP_MACRO"

# Run the analysis
echo "Starting ROOT analysis..."
root -q "$TEMP_MACRO"

# Clean up
rm -f "$TEMP_MACRO" "$TEMP_MACRO.bak"

echo ""
echo "Analysis complete!"
