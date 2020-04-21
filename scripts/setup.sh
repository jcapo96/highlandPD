if [[ -z "${HIGHLANDPATH}" ]]; then
    echo "Run the HighLAND setup script first !!!!"
    return;
fi


# and this package to the highland package hierarchy (otherwise the parameters file will not be read)
export HIGHLAND_PACKAGE_HIERARCHY=stoppingProtonAnalysis:pionAnalysis:pdBaseAnalysis:LArSoftReader:pdUtils:pdEventModel:$HIGHLAND_PACKAGE_HIERARCHY

# automatically find the highland source directory
export HIGHLANDPDPATH="$(cd -P -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
echo "HIGHLANDPDPATH=${HIGHLANDPDPATH}"

export PATH=$HIGHLANDPDPATH/bin:$PATH    
export DYLD_LIBRARY_PATH=$HIGHLANDPDPATH/lib:$DYLD_LIBRARY_PATH
export   LD_LIBRARY_PATH=$HIGHLANDPDPATH/lib:$LD_LIBRARY_PATH



#------ set the ROOT folder for every package in highland  ------
D0=$PWD
cd $HIGHLANDPDPATH/src
for D2 in *
do
    if [ ! -d $D2 ]; then
        continue
    fi
    D3=$(echo $D2 | tr "[a-z]" "[A-Z]")
    D4=${D3}ROOT
    export ${D4}=$HIGHLANDPDPATH/src/$D2
done
cd $D0

