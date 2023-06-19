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

# test it
rivet 0000.hepmc -a JET_NTUPLE_QG:MODE=DIJET:JET_R=0.4
```

## Configurations for different samples

### Z+jet (pT > 50 GeV)

```bash
rivet 0000.hepmc -a JET_NTUPLE_QG:MODE=ZJET:JET_R=0.4:ROOTFILE=ZJET.ROOT
```

### Dijet (pT > 50 GeV)

```bash
rivet 0000.hepmc -a JET_NTUPLE_QG:MODE=DIJET:JET_R=0.4:ROOTFILE=DIJET.ROOT
```

### W(->qq)Z(->vv) (pT > 500 GeV)

```bash
rivet 0000.hepmc -a JET_NTUPLE_QG:MODE=WZ:JET_R=0.8:ROOTFILE=WZ.ROOT
```

### Dijet (pT > 500 GeV)

```bash
rivet 0000.hepmc -a JET_NTUPLE_QG:MODE=DIJET:JET_R=0.8:ROOTFILE=DIJET_HPT.ROOT
```
