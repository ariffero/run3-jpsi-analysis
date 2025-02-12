This repository contains the code to fit the pT distribution of dimuons at forward rapity around the jPsi mass. It uses Run 2 Mc and Pb--Pb data at fwd rapidity.

## Code structure:
- `run-pt.py` is the interface to run the macro(s) to perform the pT fit. It can be configured using a `.yaml` file. Run it using:
  ```
  python3 run-pt.py config.yaml
  ```
- `config.yaml` is an example of the configuration file. It allows to set different parameters (see below).
- `ptFit.C` is the macro that actually perform the template fit to the pT distribution.
- `fitJPsiInPtBins.c` performs the fit to the invariant mass distribution of dimuon to find the number of jpsi in a given pT range. It is called by the script when the template for the continuum dimuons from photon-photn fusion is not used. It produces a `.txt` file with the number of jpsi in each pT range. It also saves the result of the fits in a folder called as the config file. To run it needs a "home-made library" called `savedVarInMassFits.h`.
`producePtFitHisto.C` is the macro that reads the `.txt` produced by `fitJPsiInPtBins.c` and produces the histogram with the number of jpsi as a function of pT. The histogram is used by `ptFit.C` if the template for muons from yy is not used.

## Parameters of the config file:
- `ptFitBinning` is a list that contains the edges of the bins used for the pT fit. The bins can have different width. If left empty the binning is constructed from the parameters `lowPt`, `upPt` and `nPtBins`, and the bins will all have the same width. Default value: empty list.
- `nPtBins` number of bins in pT. If `ptFitBinning` is not empty the value in the config is ignored and the values is set to the one obtained from `ptFitBinning`. Default value: 100.
- `lowPt` lower value of pT. If `ptFitBinning` is not empty the value in the config is ignored and the values is set to the one obtained from `ptFitBinning`. Default value: 0.
-  `upPt` upper value of pT. If `ptFitBinning` is not empty the value in the config is ignored and the values is set to the one obtained from `ptFitBinning`. Default value 2.
- `useContinuumTemplate` flag to decide if use a MC template for the yy->mumu process or if fit the number of jpsi vs pT instead of number of events vs pT. Flag to `true` to use the MC template. Default value: `false`.
- `binnedFit` flag to chose to do a binned fit. Feature not implement, independently of the value of this flag the fit will be a binned fit.
- `chi2Fit` fit the pT distribution using the chi2 instead of the likelihood. Set to `true` to perform a chi2 fit. Default value: `false`.
- `rapidityRange` list that contains lower and upper value of the dimuon pair rapidity. Default value: [-4, -2.5]
- `massRange_ptFit` list that contains the lower and upper value for the mass used in the pT fit. Default value: [2.85, 3.35]
- `massRange_massFit` list that contains the lower and upper value for the mass used in the mass fits. Default value: [2, 4]