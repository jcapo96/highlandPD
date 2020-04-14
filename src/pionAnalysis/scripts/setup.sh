# and this package to the highland package hierarchy (otherwise the parameters file will not be read)
export HIGHLAND_PACKAGE_HIERARCHY=pionAnalysis:$HIGHLAND_PACKAGE_HIERARCHY

# automatically find the highland source directory
export PIONANALYSISROOT="$(cd -P -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
echo "PIONANALYSISROOT=${PIONANALYSISROOT}"


export PATH=$PIONANALYSISROOT/bin:$PATH    
export DYLD_LIBRARY_PATH=$PIONANALYSISROOT/lib:$DYLD_LIBRARY_PATH
export   LD_LIBRARY_PATH=$PIONANALYSISROOT/lib:$LD_LIBRARY_PATH
