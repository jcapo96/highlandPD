#!/bin/bash

# Script to run PreliminaryK0Selection.C
# Usage: ./run_preliminary_k0_selection.sh [filename]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_FILE="$SCRIPT_DIR/PreliminaryK0Selection.C"

# Default file
DEFAULT_FILE="/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root"

# Parse command line arguments
if [[ $# -eq 0 ]]; then
    FILE="$DEFAULT_FILE"
elif [[ $# -eq 1 ]]; then
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        echo "Usage: $0 [filename]"
        echo "  filename    Input ROOT file (default: $DEFAULT_FILE)"
        echo "  -h, --help  Show this help message"
        echo ""
        echo "This script displays events with beam particle daughters that meet the selection criteria:"
        echo "- Two daughters of the beam particle whose origin is 10-50cm away from the end of the beam particle"
        echo "- These two daughters start positions are closer than 5cm"
        exit 0
    else
        FILE="$1"
    fi
else
    echo "Usage: $0 [filename]"
    echo "Use -h or --help for more information"
    exit 1
fi

# Create temporary script that loads libraries and calls the function
TEMP_SCRIPT=$(mktemp)
cat > "$TEMP_SCRIPT" << EOF
gSystem->Load("libhighland");
gSystem->Load("libhighlandPD");
.L $SCRIPT_FILE
PreliminaryK0Selection("$FILE");
EOF

echo "Running PreliminaryK0Selection with file: $FILE"
echo "Looking for events with beam particle daughters meeting the selection criteria..."
echo ""

# Execute the command
root -q "$TEMP_SCRIPT"

# Clean up temporary file
rm -f "$TEMP_SCRIPT"
