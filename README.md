# lh23-jss-rivet

## Setup

```bash
# setup LCG environment: change the arch if not using `CentOS 7`
source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/setup.sh

# build analysis
rivet-build RivetExampleAnalyses.so JET_NTUPLE_QG.cc `root-config --cflags --libs`
export RIVET_ANALYSIS_PATH=$PWD

# check the build is successful
rivet --show-analysis JET_NTUPLE_QG

# run it
rivet 0000.hepmc -a JET_NTUPLE_QG:MODE=DIJET:JET_R=0.4
```
