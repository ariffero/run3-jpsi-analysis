#include <iostream>
#include <fstream>
#include <vector>
#include "TH1D.h"
#include "TCanvas.h"
#include "string"

void producePtFitHisto(string configName = "jpsi") {

  string filename = "jpsi-" + configName + ".txt";
  std::ifstream file(filename.c_str());
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  std::vector<double> binEdges;
  std::vector<double> values;
  std::vector<double> errors;
  double val, err, xlow, xhigh;

  while (file >> val >> err >> xlow >> xhigh) {
    if (binEdges.empty() || binEdges.back() != xlow) {
      binEdges.push_back(xlow);
    }
    binEdges.push_back(xhigh);
    values.push_back(val);
    errors.push_back(err);
  }
  file.close();

  // Remove duplicate edges
  binEdges.erase(std::unique(binEdges.begin(), binEdges.end()), binEdges.end());

  int nbins = values.size();
  TH1D* hPtFit = new TH1D("hPtFit", "pt Fit", nbins, binEdges.data());
  hPtFit->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hPtFit->GetYaxis()->SetTitle("# coherent j/#psi");

  for (int i = 0; i < nbins; ++i) {
    hPtFit->SetBinContent(i + 1, values[i]);
    hPtFit->SetBinError(i + 1, errors[i]);
  }

  TCanvas* c1 = new TCanvas("c1", "Histogram", 800, 600);
  hPtFit->Draw("histo");
  hPtFit->SaveAs(Form("ptFitData-%s.root",configName.c_str()));
}