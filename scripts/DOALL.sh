# automatically find the highland source directory
HIGHLANDPATH="$(cd -P -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
echo $HIGHLANDPATH
cd $HIGHLANDPATH/highland

source scripts/cleanup.sh
source scripts/INSTALL.sh

cd ../highlandPD

source scripts/cleanup.sh
source scripts/INSTALL.sh
cd build
