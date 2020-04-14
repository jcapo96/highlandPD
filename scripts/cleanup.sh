# check that we are in the correct folder
unset found
if [ -f .git/config ]; then
    if grep highlandPD .git/config  | grep -q 'highlandPD'; then
        found=1
    fi
fi

if [ -z "${found}" ]; then
   echo "you are not in the correct folder !!! "
   return
fi

#remove everything that was created during compilation
rm -fr bin
rm -fr lib
rm -fr include
rm -fr build
rm -f HIGHLAND_PDConfig.cmake
