import yaml
import ROOT
import array
import sys
import os
import numpy as np
from pathlib import Path
import shutil
import argparse

# default values
default_values = {
  'ptFitBinning' : [],
  'lowPt': 0,
  'upPt': 2,
  'nPtBins': 100,
  'useContinuumTemplate': False,
  'binnedFit': True,
  'chi2Fit': False,
  'rapidityRange': [-4, -2.5],
  'massRange_ptFit': [2.85, 3.35],
  'massRange_massFit': [2, 4]
}

# parse arguments and read the yaml file
parser = argparse.ArgumentParser(description="Script that perform fit to the pT distribution")
parser.add_argument('yaml_file', type=str, help='Path to the YAML configuration file')

args = parser.parse_args()
yaml_file = args.yaml_file
try:
  with open(yaml_file, 'r') as file:
    data = yaml.safe_load(file)

  # take the values from the yaml config file
  ptFitBinning = data.get("ptFitBinning", default_values["ptFitBinning"])
  useContinuumTemplate = data.get("useContinuumTemplate", default_values["useContinuumTemplate"])
  binnedFit = data.get("binnedFit", default_values["binnedFit"])
  chi2Fit = data.get("chi2Fit", default_values["chi2Fit"])

  lowPt = data.get("lowPt", default_values["lowPt"])
  upPt = data.get("upPt", default_values["upPt"])
  nPtBins = data.get("nPtBins", default_values["nPtBins"])

  # check if ptFitBinning is empty and if it is fill it from
  # pT limits and number of bis
  if len(ptFitBinning)==0:
    ptFitBinning = list(np.linspace(lowPt, upPt, nPtBins + 1))
  else:
    # if ptFitBinning is not empty fix the pT limits and
    # binning to the values from ptFitBinning
    lowPt = ptFitBinning[0]
    upPt = ptFitBinning[len(ptFitBinning)-1]
    nPtBins = len(ptFitBinning)-1

  # Convert Python list to ROOT std::vector
  ptBinning = ROOT.std.vector("double")(ptFitBinning)

  rapidityRange = data.get("rapidityRange", default_values["rapidityRange"])
  rapidity = array.array("d", rapidityRange)

  massRange_ptFit = data.get("massRange_ptFit", default_values["massRange_ptFit"])
  mass_ptFit = array.array("d", massRange_ptFit)

  massRange_massFit = data.get("massRange_massFit", default_values["massRange_massFit"])
  mass_massFit = array.array("d", massRange_massFit)

  configName = yaml_file.strip('.yaml')
  # if the pT fit requires a histo with jpsi from mass fits
  # check if the histogram corresponding a certain config already exists
  # and has not been modified after modifying the config
  if not useContinuumTemplate:
    jpsiHisto = 'ptFitData-' + configName + '.root'
    jpsiTxt = 'jpsi-' + configName + '.txt'

    time_yaml = Path(yaml_file).stat().st_mtime
    
    time_histo = 0
    if os.path.exists(jpsiHisto):
      time_histo = Path(jpsiHisto).stat().st_mtime
    
    if os.path.exists(jpsiHisto) and time_histo > time_yaml:
      print("The histo already exists, not doing mass fits.")
    else:
      if os.path.exists(jpsiTxt):
        os.remove(jpsiTxt)
      if Path(configName).exists() and Path(configName).is_dir():
        shutil.rmtree(Path(configName)) 
      # do the mass fits
      for ptId in range(1, len(ptBinning)):
        ROOT.gROOT.LoadMacro("fitJPsiInPtBins.c+")
        ROOT.fitJPsiInPtBins(ptBinning, ptId, mass_massFit, rapidity,configName)
        # unload the compiled macro to avoid issue of variables with the same name
        ROOT.gROOT.ProcessLine(".U fitJPsiInPtBins_c.so")

      # produce the histo with jpsi vs pT
      ROOT.gROOT.LoadMacro("producePtFitHisto.c+")
      ROOT.producePtFitHisto(configName)
      ROOT.gROOT.ProcessLine(".U producePtFitHisto_c.so")

  # do the pT fit
  ROOT.gROOT.LoadMacro("ptFit.C+")
  ROOT.ptFit(ptBinning, useContinuumTemplate, binnedFit, chi2Fit, mass_ptFit, rapidity, configName)

  # keep ROOT open after running
  input("\nPress Enter to exit...")

except FileNotFoundError:
  print(f"Error: The file '{args.yaml_file}' was not found.")
except yaml.YAMLError as e:
  print(f"Error parsing YAML file: {e}")
except Exception as e:
  print(f"Error loading configuration file: {e}")