if [[ -z "${HIGHLANDPATH}" ]]; then
    echo "Run the HighLAND setup script first !!!!"
    return
fi

# and this package to the highland package hierarchy (otherwise the parameters file will not be read)
#export HIGHLAND_PACKAGE_HIERARCHY=secondaryKaonAnalysis:pionAnalysis:stoppingProtonAnalysis:pdStudies:pdSystematics:pdBaseAnalysis:LArSoftReader:pdUtils:pdEventModel:$HIGHLAND_PACKAGE_HIERARCHY

# Determine this script's directory robustly for both bash and zsh, then set HIGHLANDPDPATH
__hlpd_script_path=""
if [[ -n "${BASH_SOURCE[0]}" ]]; then
    __hlpd_script_path="${BASH_SOURCE[0]}"
elif [[ -n "${ZSH_VERSION}" ]]; then
    # zsh: %N expands to the path of the sourced script
    __hlpd_script_path="${(%):-%N}"
else
    __hlpd_script_path="$0"
fi
__hlpd_script_dir="$(cd -P -- "$(dirname -- "${__hlpd_script_path}")" 2>/dev/null && pwd -P)"
if [[ -z "${__hlpd_script_dir}" ]]; then
    echo "Error: could not determine highlandPD setup script directory"
    return 1
fi
export HIGHLANDPDPATH="$(cd -P -- "${__hlpd_script_dir}/.." && pwd -P)"
echo "HIGHLANDPDPATH=${HIGHLANDPDPATH}"

export PATH="${HIGHLANDPDPATH}/bin:${PATH}"
export DYLD_LIBRARY_PATH="${HIGHLANDPDPATH}/lib:${DYLD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${HIGHLANDPDPATH}/lib:${LD_LIBRARY_PATH}"

export CMAKE_PREFIX_PATH="${HIGHLANDPDPATH}/cmake/modules:${CMAKE_PREFIX_PATH}"

#------ set the ROOT folder for every package in highlandPD ------
D0="$PWD"
if [[ -d "${HIGHLANDPDPATH}/src" ]]; then
    cd "${HIGHLANDPDPATH}/src" || return 1
    for D2 in *; do
        if [[ ! -d "$D2" ]]; then
            continue
        fi
        D3=$(echo "$D2" | tr "[a-z]" "[A-Z]")
        D4=${D3}ROOT
        export ${D4}="${HIGHLANDPDPATH}/src/${D2}"
    done
    cd "$D0" || return 1
else
    echo "Warning: ${HIGHLANDPDPATH}/src not found; package ROOT variables not set"
fi

