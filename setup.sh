cd /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_5_0
cmsenv
cd -;
export LANG=en_US.UTF-8

[ -d myenv ] || virtualenv myenv
source myenv/bin/activate

[[ -d CORE ]] || {
    git clone git@github.com:cmstas/CORE.git
    cd CORE
    make -j10
    cd -
}

# pip install --user uproot
