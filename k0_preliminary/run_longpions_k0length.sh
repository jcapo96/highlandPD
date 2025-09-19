#!/bin/bash
# run_longpions_k0length.sh
# Wrapper script to run LongPionsK0Length.C

# Default file
DEFAULT_FILE="/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root"

# Use provided file or default
FILE=${1:-$DEFAULT_FILE}

echo "Running LongPionsK0Length.C with file: $FILE"
echo "This script calculates distances between K0 parent end position and π+/- start positions"
echo "Uses truth information to find K0 and parent, but calculates distances using RECO information"
echo "Only processes events with π+ and π- daughters that travel at least 10cm each"
echo ""

# Run the ROOT macro
root -q -b "LongPionsK0Length.C(\"$FILE\")"
