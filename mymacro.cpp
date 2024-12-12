#include <TApplication.h>
#include <TBrowser.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TStyle.h>

#include "particle.hpp"

//
Double_t Gauss(Double_t *x, Double_t *par) {
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. /
                                     par[2] / par[2]);
  return val;
}

void mymacro() {
  Int_t BINS = 1E3;

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
  TH1F *hDiff1 = new TH1F("hDiff1", "Diff1", BINS, 0, 7);
  hDiff1->Sumw2();
  TH1F *hDiff2 = new TH1F("hDiff2", "Diff2", BINS, 0, 7);
  hDiff2->Sumw2();

  int NBins;

  double expectedProportions[] = {0.4, 0.4, 0.05, 0.05, 0.045, 0.045, 0.01};
  int nTypes = 7;

  int Entries;
  Entries = hParticleTypes->GetEntries();
  if (hParticleTypes && Entries == 1E7)
    std::cout << "Ok, hParticle has " << Entries << " entries; " << std::endl;
  else
    std::cout << "Something with hParticle is wrong " << std::endl;

  NBins = hParticleTypes->GetNbinsX();
  for (int bin = 1; bin <= NBins; ++bin) {
    double content = hParticleTypes->GetBinContent(bin);
    double error = hParticleTypes->GetBinError(bin);
    double expectedContent = expectedProportions[bin - 1] * Entries;
    double z = (content - expectedContent) / error;
    std::cout << "Bin[" << bin << "] has " << content << " +/- " << error
              << " entries; the standard score z:\nz = " << z << std::endl;
  }

  double chisquare;
  double NDF;
  double probFit;

  Entries = hPolarAngles->GetEntries();
  NBins = hPolarAngles->GetNbinsX();
  if (hPolarAngles && hPolarAngles->GetEntries() == 1E7)
    std::cout << "Ok, hPolarAngles has " << hPolarAngles->GetEntries()
              << " entries; " << std::endl;
  else
    std::cout << "Something with hPolarAngles is wrong " << std::endl;

  TF1 *fitPolar = new TF1("Fit Polar Angles", "[0]", 0, TMath::Pi());
  fitPolar->SetParameter(0, Entries / NBins);
  hPolarAngles->Fit(fitPolar, "APE");
  chisquare = fitPolar->GetChisquare();
  NDF = fitPolar->GetNDF();
  probFit = fitPolar->GetProb();
  std::cout << "chiSquare/NDf= " << chisquare / NDF << std::endl;
  if (probFit > 0.05)
    std::cout << "The observed distribution is consistent with the expected "
                 "one.\nThe probability to get a higher chi square than "
              << chisquare << " is: " << probFit << std::endl;
  else
    std::cout << "The observed distribution is inconsistent with the expected "
                 "one.\nThe probability to get a higher chi square than "
              << chisquare << " is: " << probFit << std::endl;

  Entries = hAzimuthalAngles->GetEntries();
  NBins = hAzimuthalAngles->GetNbinsX();
  if (hAzimuthalAngles && hAzimuthalAngles->GetEntries() == 1E7)
    std::cout << "Ok, hAzimuthalAngles has " << hAzimuthalAngles->GetEntries()
              << " entries; " << std::endl;
  else
    std::cout << "Something with hAzimuthalAngles is wrong " << std::endl;

  TF1 *fitAzhimutal =
      new TF1("Fit Azhimutal Angles", "[0]", 0, 2 * TMath::Pi());
  fitAzhimutal->SetParameter(0, Entries / NBins);
  hAzimuthalAngles->Fit(fitAzhimutal, "APE");
  chisquare = fitPolar->GetChisquare();
  NDF = fitAzhimutal->GetNDF();
  probFit = fitAzhimutal->GetProb();
  std::cout << "Fit parameter: " << fitAzhimutal->GetParameter(0) << std::endl;
  std::cout << "chiSquare/NDf = " << chisquare / NDF << std::endl;
  if (probFit > 0.05)
    std::cout << "The observed distribution is consistent with the expected "
                 "one.\nThe probability to get a higher chi square than "
              << chisquare << " is: " << probFit << std::endl;
  else
    std::cout << "The observed distribution is inconsistent with the expected "
                 "one.\nThe probability to get a higher chi square than "
              << chisquare << " is: " << probFit << std::endl;

  Entries = hImpulse->GetEntries();
  if (hImpulse && hImpulse->GetEntries() == 1E7)
    std::cout << "Ok, hImpulse has " << hImpulse->GetEntries() << " entries; "
              << std::endl;
  else
    std::cout << "Something with hImpulse is wrong " << std::endl;

  TF1 *fitImpulse = new TF1("Fit Impulse", "expo", 0, 10);
  hImpulse->Fit(fitImpulse, "APE");
  chisquare = fitImpulse->GetChisquare();
  NDF = fitImpulse->GetNDF();
  probFit = fitImpulse->GetProb();
  std::cout << "Fit parameter: [0] = " << fitImpulse->GetParameter(0)
            << "; [1] = " << fitImpulse->GetParameter(1) << std::endl;
  std::cout << "chiSquare/NDf = " << chisquare / NDF << std::endl;
  if (probFit > 0.05)
    std::cout << "The observed distribution is consistent with the expected "
                 "one.\nThe probability to get a higher chi square than "
              << chisquare << " is: " << probFit << std::endl;
  else
    std::cout << "The observed distribution is inconsistent with the expected "
                 "one.\nThe probability to get a higher chi square than "
              << chisquare << " is: " << probFit << std::endl;

  Entries = hTransverseImpulse->GetEntries();
  if (hTransverseImpulse && hTransverseImpulse->GetEntries() == 1E7)
    std::cout << "Ok, hTransverseImpulse has "
              << hTransverseImpulse->GetEntries() << " entries; " << std::endl;
  else
    std::cout << "Something with hTransverseImpulse is wrong " << std::endl;

  Entries = hEnergy->GetEntries();
  if (hEnergy && hEnergy->GetEntries() == 1E7)
    std::cout << "Ok, hEnergy has " << hTransverseImpulse->GetEntries()
              << " entries; " << std::endl;
  else
    std::cout << "Something with hEnergy is wrong " << std::endl;

  Entries = hInvMassAll->GetEntries();
  if (hInvMassAll)
    std::cout << "Ok, hInvMassAll has " << hInvMassAll->GetEntries()
              << " entries; " << std::endl;
  else
    std::cout << "Something with hInvMassAll is wrong " << std::endl;
  if (hInvMassCon)
    std::cout << "Ok, hInvMassCon has " << hInvMassCon->GetEntries()
              << " entries; " << std::endl;
  else
    std::cout << "Something with hInvMassCon is wrong " << std::endl;
  if (hInvMassDis)
    std::cout << "Ok, hInvMassDis has " << hInvMassDis->GetEntries()
              << " entries; " << std::endl;
  else
    std::cout << "Something with hInvMassDis is wrong " << std::endl;
  if (hInvMassConPiK)
    std::cout << "Ok, hInvMassConPiK has " << hInvMassConPiK->GetEntries()
              << " entries; " << std::endl;
  else
    std::cout << "Something with hInvMassConPiK is wrong " << std::endl;
  if (hInvMassDisPiK)
    std::cout << "Ok, hInvMassDisPiK has " << hInvMassDisPiK->GetEntries()
              << " entries; " << std::endl;
  else
    std::cout << "Something with hInvMassDisPiK is wrong " << std::endl;
  if (hBenchmark)
    std::cout << "Ok, hBenchmark has " << hBenchmark->GetEntries()
              << " entries; " << std::endl;
  else
    std::cout << "Something with hBenchmark is wrong " << std::endl;
  /*if (hKStar)
    std::cout << "Ok, hKStar has " << hKStar->GetEntries() << " entries; "
              << std::endl;
  else
    std::cout << "Something with hKStar is wrong " << std::endl;*/

  hDiff1->Add(hInvMassDis, hInvMassCon, 1, -1);
  hDiff2->Add(hInvMassDisPiK, hInvMassConPiK, 1, -1);

  TF1 *fitKStar1 = new TF1("fitKStar1", Gauss, 0, 7, 3);
  fitKStar1->SetNpx(2000);
  fitKStar1->SetParameters(hDiff1->GetMaximum(), 0.89166, 0.05);
  hDiff1->Fit(fitKStar1);
  double extractedMass1 = fitKStar1->GetParameter(1);
  double MasseError1 = fitKStar1->GetParError(1);
  double extractedWidth1 = fitKStar1->GetParameter(2);
  double WidthError1 = fitKStar1->GetParError(2);

  std::cout << "Massa del K* estratta da hDiff1: " << extractedMass1 << " +/- "
            << MasseError1 << std::endl;
  std::cout << "Larghezza del K* estratta da hDiff1: " << extractedWidth1
            << " +/- " << WidthError1 << std::endl;
  std::cout << "Chi square/ NDF of hDiff1: "
            << fitKStar1->GetChisquare() / fitKStar1->GetNDF() << std::endl;
  std::cout << "the probability"
            << std::endl;  // vediamo da come abbiamo fatto prima
  std::cout << "The maximum content is: " << hDiff1->GetMaximum() << std::endl;

  TF1 *fitKStar2 = new TF1("fitKStar2", Gauss, 0, 7, 3);
  fitKStar2->SetNpx(2000);
  fitKStar2->SetParameters(hDiff2->GetMaximum(), 0.89166, 0.05);
  hDiff2->Fit(fitKStar2, "S");
  double extractedMass2 = fitKStar2->GetParameter(1);
  double extractedWidth2 = fitKStar2->GetParameter(2);

  std::cout << "Massa del K* estratta da hDiff2: " << extractedMass2
            << std::endl;
  std::cout << "Larghezza del K* estratta da hDiff2: " << extractedWidth2
            << std::endl;
  std::cout << "Chi square/ NDF of hDiff1: "
            << fitKStar2->GetChisquare() / fitKStar2->GetNDF() << std::endl;
  std::cout << "the probability";  // vediamo da come abbiamo fatto prima

  TF1 *fitBenchmark = new TF1("fitBenchmark", Gauss, 0.7, 1.1, 3);
  fitBenchmark->SetParameters(hBenchmark->GetMaximum(), 0.89166, 0.05);
  fitBenchmark->SetNpx(2000);
  hBenchmark->Fit(fitBenchmark, "S, L");

  // zona cosmetica

  hParticleTypes->SetName("Particle Types");
  hParticleTypes->GetXaxis()->SetTitle("Particle Types");
  hParticleTypes->GetYaxis()->SetTitle("Counts");
  hParticleTypes->SetFillColor(kBlue);
  hParticleTypes->SetLineWidth(0);
  hParticleTypes->SetMarkerSize(0);
  hParticleTypes->GetXaxis()->SetBinLabel(1, "Pi+");
  hParticleTypes->GetXaxis()->SetBinLabel(2, "Pi-");
  hParticleTypes->GetXaxis()->SetBinLabel(3, "K+");
  hParticleTypes->GetXaxis()->SetBinLabel(4, "K-");
  hParticleTypes->GetXaxis()->SetBinLabel(5, "P+");
  hParticleTypes->GetXaxis()->SetBinLabel(6, "P-");
  hParticleTypes->GetXaxis()->SetBinLabel(7, "K*");

  hPolarAngles->SetName("Polar angles");
  hPolarAngles->GetXaxis()->SetTitle("Polar Angle (rad)");
  hPolarAngles->GetYaxis()->SetTitle("Counts");
  hPolarAngles->SetLineColor(kBlue);
  hPolarAngles->SetLineWidth(2);
  hPolarAngles->GetYaxis()->SetRangeUser(9000, 11000);

  hAzimuthalAngles->SetName("Azimutal angles");
  hAzimuthalAngles->GetXaxis()->SetTitle("Azimuthal Angle (rad)");
  hAzimuthalAngles->GetYaxis()->SetTitle("Counts");
  hAzimuthalAngles->SetLineColor(kBlue);
  hAzimuthalAngles->SetLineWidth(2);
  hAzimuthalAngles->GetYaxis()->SetRangeUser(9000, 11000);

  hImpulse->SetName("Impulse");
  hImpulse->GetXaxis()->SetTitle("Impulse (GeV/c)");
  hImpulse->GetYaxis()->SetTitle("Counts");
  hImpulse->SetLineColor(kBlue);
  hImpulse->SetLineWidth(2);

  hTransverseImpulse->SetName("Trasverse impulse");
  hTransverseImpulse->GetXaxis()->SetTitle("Transverse Impulse (GeV/c)");
  hTransverseImpulse->GetYaxis()->SetTitle("Counts");
  hTransverseImpulse->SetLineColor(kBlue);
  hTransverseImpulse->SetLineWidth(2);

  hEnergy->SetName("Energy");
  hEnergy->GetXaxis()->SetTitle("Energy (GeV)");
  hEnergy->GetYaxis()->SetTitle("Counts");
  hEnergy->SetLineColor(kBlue);
  hEnergy->SetLineWidth(2);

  hInvMassAll->SetName("InvMass between all particles");
  hInvMassAll->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hInvMassAll->GetYaxis()->SetTitle("Counts");
  hInvMassAll->SetLineColor(kBlue);
  hInvMassAll->SetLineWidth(2);

  hInvMassCon->SetName("InvMass between concordant charge particles");
  hInvMassCon->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hInvMassCon->GetYaxis()->SetTitle("Counts");
  hInvMassCon->SetLineColor(kBlue);
  hInvMassCon->SetLineWidth(2);

  hInvMassDis->SetName("InvMass between discordant charge particles");
  hInvMassDis->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hInvMassDis->GetYaxis()->SetTitle("Counts");
  hInvMassDis->SetLineColor(kBlue);
  hInvMassDis->SetLineWidth(2);

  hInvMassConPiK->SetName(
      "InvMass between concordant charge Pi and K particles");
  hInvMassConPiK->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hInvMassConPiK->GetYaxis()->SetTitle("Counts");
  hInvMassConPiK->SetLineColor(kBlue);
  hInvMassConPiK->SetLineWidth(2);

  hInvMassDisPiK->SetName(
      "InvMass between discordant charge Pi and K particles");
  hInvMassDisPiK->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hInvMassDisPiK->GetYaxis()->SetTitle("Counts");
  hInvMassDisPiK->SetLineColor(kBlue);
  hInvMassDisPiK->SetLineWidth(2);

  hBenchmark->SetName("Benchmark");
  hBenchmark->GetXaxis()->SetTitle("Benchmark Variable");
  hBenchmark->GetYaxis()->SetTitle("Counts");
  hBenchmark->SetLineColor(kBlue);
  hBenchmark->SetLineWidth(2);
  hBenchmark->GetXaxis()->SetRangeUser(0.6, 1.2);

  // Custom range for Diff histograms based on observed data ranges
  hDiff1->SetName("Diff1");
  hDiff1->GetXaxis()->SetTitle("Invariant Mass Difference (GeV/c^2)");
  hDiff1->GetYaxis()->SetTitle("Counts");
  hDiff1->SetLineColor(kBlue);
  hDiff1->SetLineWidth(2);
  hDiff1->GetXaxis()->SetRangeUser(0.6, 1.2);
  fitKStar1->SetLineColor(kRed);
  fitKStar1->SetLineWidth(4);

  hDiff2->SetName("Diff2");
  hDiff2->GetXaxis()->SetTitle("Invariant Mass Difference (GeV/c^2)");
  hDiff2->GetYaxis()->SetTitle("Counts");
  hDiff2->SetLineColor(kBlue);
  hDiff2->SetLineWidth(2);
  hDiff2->GetXaxis()->SetRangeUser(0.6, 1.2);

  TCanvas *c1 =
      new TCanvas("c1", "Particle and its cinematic properties", 750, 750);
  c1->Divide(2, 2);
  c1->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptStat(11);
  hParticleTypes->Draw("BAR2");
  c1->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptFit(1111);
  hImpulse->Draw("HIST");
  fitImpulse->Draw("SAME");
  c1->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  hPolarAngles->Draw("HIST");
  fitPolar->Draw("SAME");
  c1->cd(4);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  hAzimuthalAngles->Draw("HIST");
  fitAzhimutal->Draw("SAME");
  c1->SaveAs("Particles_Properties.pdf");
  c1->SaveAs("Particles_Properties.C");
  c1->SaveAs("Particles_Properties.root");

  TCanvas *c2 = new TCanvas("c2", "Invariant mass difference", 500, 800);
  c2->Divide(1, 3);
  c2->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  hBenchmark->Draw("HIST");
  fitBenchmark->Draw("SAME");
  c2->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  hDiff1->Draw("HIST");
  fitKStar1->Draw("SAME");
  c2->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  hDiff2->Draw("HIST");
  fitKStar2->Draw("SAME");
  c2->SaveAs("Invariant_mass_difference.pdf");
  c2->SaveAs("Invariant_mass_difference.C");
  c2->SaveAs("Invariant_mass_difference.root");

  TCanvas *c3 = new TCanvas("c3", "TransverseImpulse", 375, 375);
  gStyle->SetOptFit(0000);
  hTransverseImpulse->Draw("HIST");
  c3->SaveAs("TransverseImpulse.pdf");
  c3->SaveAs("TransverseImpulse.C");
  c3->SaveAs("TransverseImpulse.root");

  TCanvas *c4 = new TCanvas("c4", "Energy", 375, 375);
  hEnergy->Draw("HIST");
  c4->SaveAs("Energy.pdf");
  c4->SaveAs("Energy.C");
  c4->SaveAs("Energy.root");

  TCanvas *c5 = new TCanvas("c5", "InvMassAll", 375, 375);
  hInvMassAll->Draw("HIST");
  c5->SaveAs("InvMassAll.pdf");
  c5->SaveAs("InvMassAll.C");
  c5->SaveAs("InvMassAll.root");

  TCanvas *c6 = new TCanvas("c6", "Invariant Mass", 750, 750);
  c6->Divide(2, 2);
  c6->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  hInvMassCon->Draw("HIST");
  c6->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  hInvMassDis->Draw("HIST");
  c6->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  hInvMassConPiK->Draw("HIST");
  c6->cd(4);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  hInvMassDisPiK->Draw("HIST");
  c6->SaveAs("InvMass.pdf");
  c6->SaveAs("InvMass.C");
  c6->SaveAs("InvMass.root");
}
