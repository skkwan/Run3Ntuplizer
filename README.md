# Run3Ntuplizer

Based on [Run3Ntuplizer](https://github.com/isobelojalvo/Run3Ntuplizer/tree/2020_Mar_5). Emulate Run 3 and create
n-tuples of L1 taus.

## Installation

```
cmsrel CMSSW_10_6_0_pre4
cd CMSSW_10_6_0_pre4/src
cmsenv
git cms-init
git cms-merge-topic isobelojalvo:run3-dev-$CMSSW_VERSION
cd L1Trigger
git clone git@github.com:isobelojalvo/L1TCaloSummary.git
git clone git@github.com:skkwan/Run3Ntuplizer.git
cd ..
scram b -j 8

cd L1Trigger/Run3Ntuplizer/test
cmsRun test-Analyzer.py
```