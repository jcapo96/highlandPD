#!/bin/bash

# Script to run CollectK0Masses.C
# Usage: ./run_collect_masses.sh [filename]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_FILE="$SCRIPT_DIR/CollectK0Masses.C"

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
CollectK0Masses("$FILE");
EOF

echo "Running CollectK0Masses with file: $FILE"
echo ""

# Execute the command
root -q "$TEMP_SCRIPT"

# Clean up temporary file
rm -f "$TEMP_SCRIPT"
