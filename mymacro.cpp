#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom.h>

#include "particle.hpp"

void mymacro() {
  TFile *file = new TFile("histograms.root");
  TH1F *hParticleTypes = (TH1F *)file->Get("h1");
  TH1F *hPolarAngles = (TH1F *)file->Get("h2");
  TH1F *hAzimuthalAngles = (TH1F *)file->Get("h3");
  TH1F *hImpulse = (TH1F *)file->Get("h4");
  TH1F *hTransverseImpulse = (TH1F *)file->Get("h5");
  TH1F *hEnergy = (TH1F *)file->Get("h6");
  TH1F *hInvMassAll = (TH1F *)file->Get("h7");
  TH1F *hInvMassCon = (TH1F *)file->Get("h8");
  TH1F *hInvMassDis = (TH1F *)file->Get("h9");
  TH1F *hInvMassConPiK = (TH1F *)file->Get("h10");
  TH1F *hInvMassDisPiK = (TH1F *)file->Get("h11");
  TH1F *hBenchmark = (TH1F *)file->Get("h12");
  TH1F *hKStar = (TH1F *)file->Get("h13");
  int NBins;

  double expectedProportions[] = {
      0.4,   0.4,   0.05, 0.05,
      0.045, 0.045, 0.01};  // Proporzioni attese per 3 tipi di particelle
  int nTypes = 7;

  int Entries = hParticleTypes->GetEntries();
  if (hParticleTypes && Entries == 1E5)
    std::cout << "Ok, hParticle has " << Entries << " entries; " << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;

  NBins = hParticleTypes->GetNbinsX();
  for (int bin = 1; bin <= NBins; ++bin) {
    double content = hParticleTypes->GetBinContent(bin);
    double error = hParticleTypes->GetBinError(bin);
    double observedProportion = content / Entries;
    double expectedProportion = expectedProportions[bin - 1];
    double z = (observedProportion - expectedProportion) / (error / Entries); 
    std::cout<< "Bin[ " << bin << " ] has a standard normal variable: z = " << z << std::endl; 
  }

  if (hPolarAngles && hPolarAngles->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;

  NBins = hParticleTypes->GetNbinsX();
  for (int bin = 1; bin <= NBins; ++bin) {
    double content = hParticleTypes->GetBinContent(bin);
    double error = hParticleTypes->GetBinError(bin);
  }
  if (hAzimuthalAngles && hAzimuthalAngles->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
  if (hImpulse && hImpulse->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
  if (hTransverseImpulse && hTransverseImpulse->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
  if (hEnergy && hEnergy->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
  if (hInvMassAll && hInvMassAll->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
  if (hInvMassCon && hInvMassCon->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
  if (hInvMassDis && hInvMassDis->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
  if (hInvMassConPiK && hInvMassConPiK->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
  if (hInvMassDisPiK && hInvMassDisPiK->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
  if (hBenchmark && hBenchmark->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
  if (hKStar && hKStar->GetEntries() == 1E5)
    std::cout << "Ok" << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;
}