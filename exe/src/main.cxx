#include "landgausFit.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "fitFunctionClass.hh"
int main(int argc, char **argv){
  
  TApplication app("App", &argc, argv);




  TFile* f = new TFile("C:/slac_data/abc/run_0712_0935/test_12.root");
  f->ls();
  TTree* signal = (TTree*)f->Get("signal;1");
  landgausFit fit;
  
  fit.setStartAmplitude(99);
  fit.setLimits_Amplitude(90, 100);
  fit.setLimits_GaussSigma(0, 200);
  fit.setLimits_LandauMP(0, 500);
  fit.setLimits_LandauSigma(0, 200);
  fit.setStartGaussSigma(27);
  fit.setStartLandauSigma(7);
  fit.setFitRange(0,250);
  fit.setFitOptions("");
  TCanvas c;
  std::string name,cutpara;
  TFile* outf = new TFile("C:/slac_data/abc/run_0712_0935/test_12_out.root","recreate");
  TTree * tree = new TTree("lfits","landauFitParameters");

  Double_t m_landau_mean,m_landau_sigma,m_gaus_sigma,m_amplitude,chiOverNDF;
  Int_t channelNR;
  tree->Branch("landau_mean", &m_landau_mean, "landau_mean/D");
  tree->Branch("landau_sigma", &m_landau_sigma, "landau_sigma/D");
  tree->Branch("gaus_sigma", &m_gaus_sigma, "gaus_sigma/D");
  tree->Branch("amplitude", &m_amplitude, "amplitude/D");
  tree->Branch("chiOverNDF", &chiOverNDF, "chiOverNDF/D");
  tree->Branch("channelNR", &channelNR, "channelNR/I");
  for (Int_t i = 0; i < 102;++i)
  {
    channelNR = i;
    cutpara = "isolating==5&&channel==" + std::to_string(i);
    name = "effiVsThr_for_channel_" + std::to_string(i) + ".png";
    signal->Draw("effi4:threshold", cutpara.c_str(), "*");
    TGraph * g = (TGraph*) c.GetPrimitive("Graph");


    fit(g);
    fit.printResults();
    fit.DrawfitFunction();
    c.Print(name.c_str());
    m_landau_mean = fit.getLandauMostProbable();
    m_landau_sigma = fit.getLandauSigma();
    m_gaus_sigma = fit.getGaussSigma();
    m_amplitude = fit.getAmplitude();
    chiOverNDF = fit.getChiSqare() / fit.GetNDF();
    tree->Fill();
  }
  outf->Write();
  outf->Close();
 // app.Run();

}