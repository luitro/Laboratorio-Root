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

#include "particle.hpp"

// Gaussiana 1-D definita nella funzione utente
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
  TH1F *hKStar = (TH1F *)file->Get("h13");
  TH1F *hDiff1 = new TH1F("hDiff1", "Diff1", BINS, 0, 7);
  hDiff1->Sumw2();
  TH1F *hDiff2 = new TH1F("hDiff2", "Diff2", BINS, 0, 7);
  hDiff2->Sumw2();

  int NBins;

  double expectedProportions[] = {
      0.4,   0.4,   0.05, 0.05,
      0.045, 0.045, 0.01};  // Proporzioni attese per 3 tipi di particelle
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
    /*double observedProportion = content / Entries;
    double expectedProportion = expectedProportions[bin - 1];
    double z = (observedProportion - expectedProportion) / (error / Entries);
    std::cout << "Bin[ " << bin
              << " ] has a standard normal variable: z = " << z << std::endl;*/
    double expectedContent = expectedProportions[bin - 1] * Entries;
    double z = (content - expectedContent) / error;
    std::cout << "Bin[" << bin << "] has " << content
              << " entries and a standard normal variable: z = " << z
              << "\nContent: " << content << " +/- " << error << std::endl;
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
  if (hKStar)
    std::cout << "Ok, hKStar has " << hKStar->GetEntries() << " entries; "
              << std::endl;
  else
    std::cout << "Something with hKStar is wrong " << std::endl;

  hDiff1->Add(hInvMassDis, hInvMassCon, 1, -1);
  hDiff2->Add(hInvMassDisPiK, hInvMassConPiK, 1, -1);

  TF1 *fitKStar1 = new TF1("fitKStar1", Gauss, 0, 7, 3);
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
  std::cout << "the probability";  // vediamo da come abbiamo fatto prima


  TF1 *fitKStar2 = new TF1("fitKStar2", Gauss, 0, 7, 3);
  fitKStar2->SetParameters(hDiff2->GetMaximum(), 0.89166, 0.05);
  hDiff2->Fit(fitKStar2);
  double extractedMass2 = fitKStar2->GetParameter(1);
  double extractedWidth2 = fitKStar2->GetParameter(2);

  std::cout << "Massa del K* estratta da hDiff2: " << extractedMass2
            << std::endl;
  std::cout << "Larghezza del K* estratta da hDiff2: " << extractedWidth2
            << std::endl;
  std::cout << "Chi square/ NDF of hDiff1: "
            << fitKStar2->GetChisquare() / fitKStar2->GetNDF() << std::endl;
  std::cout << "the probability";  // vediamo da come abbiamo fatto prima

  // zona cosmetica

  hParticleTypes->GetXaxis()->SetTitle("Particle Types");
  hParticleTypes->GetYaxis()->SetTitle("Counts");
  hParticleTypes->SetFillColor(kBlue);
  hParticleTypes->SetLineWidth(2);

  hPolarAngles->GetXaxis()->SetTitle("Polar Angle (rad)");
  hPolarAngles->GetYaxis()->SetTitle("Counts");
  hPolarAngles->SetLineColor(kBlue);
  hPolarAngles->SetLineWidth(2);
  hPolarAngles->GetYaxis()->SetRangeUser(9000, 11000);

  hAzimuthalAngles->GetXaxis()->SetTitle("Azimuthal Angle (rad)");
  hAzimuthalAngles->GetYaxis()->SetTitle("Counts");
  hAzimuthalAngles->SetLineColor(kBlue);
  hAzimuthalAngles->SetLineWidth(2);

  hImpulse->GetXaxis()->SetTitle("Impulse (GeV/c)");
  hImpulse->GetYaxis()->SetTitle("Counts");
  hImpulse->SetLineColor(kBlue);
  hImpulse->SetLineWidth(2);

  hTransverseImpulse->GetXaxis()->SetTitle("Transverse Impulse (GeV/c)");
  hTransverseImpulse->GetYaxis()->SetTitle("Counts");
  hTransverseImpulse->SetLineColor(kBlue);
  hTransverseImpulse->SetLineWidth(2);

  hEnergy->GetXaxis()->SetTitle("Energy (GeV)");
  hEnergy->GetYaxis()->SetTitle("Counts");
  hEnergy->SetLineColor(kBlue);
  hEnergy->SetLineWidth(2);

  hInvMassAll->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hInvMassAll->GetYaxis()->SetTitle("Counts");
  hInvMassAll->SetLineColor(kBlue);
  hInvMassAll->SetLineWidth(2);

  hInvMassCon->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hInvMassCon->GetYaxis()->SetTitle("Counts");
  hInvMassCon->SetLineColor(kBlue);
  hInvMassCon->SetLineWidth(2);

  hInvMassDis->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hInvMassDis->GetYaxis()->SetTitle("Counts");
  hInvMassDis->SetLineColor(kBlue);
  hInvMassDis->SetLineWidth(0);
  hInvMassDis->SetFillColor(kBlue);

  hInvMassConPiK->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hInvMassConPiK->GetYaxis()->SetTitle("Counts");
  hInvMassConPiK->SetLineColor(kBlue);
  hInvMassConPiK->SetLineWidth(2);

  hInvMassDisPiK->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hInvMassDisPiK->GetYaxis()->SetTitle("Counts");
  hInvMassDisPiK->SetLineColor(kBlue);
  hInvMassDisPiK->SetLineWidth(2);

  hBenchmark->GetXaxis()->SetTitle("Benchmark Variable");
  hBenchmark->GetYaxis()->SetTitle("Counts");
  hBenchmark->SetLineColor(kBlue);
  hBenchmark->SetLineWidth(2);

  hKStar->GetXaxis()->SetTitle("K* Mass (GeV/c^2)");
  hKStar->GetYaxis()->SetTitle("Counts");
  hKStar->SetLineColor(kBlue);
  hKStar->SetLineWidth(2);

  // Custom range for Diff histograms based on observed data ranges
  hDiff1->GetXaxis()->SetTitle("Invariant Mass Difference (GeV/c^2)");
  hDiff1->GetYaxis()->SetTitle("Counts");
  hDiff1->SetLineColor(kBlue);
  hDiff1->SetLineWidth(2);

  hDiff2->GetXaxis()->SetTitle("Invariant Mass Difference (GeV/c^2)");
  hDiff2->GetYaxis()->SetTitle("Counts");
  hDiff2->SetLineColor(kBlue);
  hDiff2->SetLineWidth(2);

  // Set vertical scale (if needed for clarity)
  // hDiff1->GetYaxis()->SetRangeUser(hDiff1->GetMinimum() * 1.2,
  // hDiff1->GetMaximum() * 1.2);
  // hDiff2->GetYaxis()->SetRangeUser(hDiff2->GetMinimum() * 1.2,
  // hDiff2->GetMaximum() * 1.2);

  TCanvas *c1 = new TCanvas("c1", "Canvas for hParticleTypes", 800, 600);
  hParticleTypes->Draw();
  c1->SaveAs("hParticleTypes.pdf");
  c1->SaveAs("hParticleTypes.C");
  c1->SaveAs("hParticleTypes.root");

  TCanvas *c2 = new TCanvas("c2", "Canvas for hPolarAngles", 800, 600);
  hPolarAngles->Draw();
  c2->SaveAs("hPolarAngles.pdf");
  c2->SaveAs("hPolarAngles.C");
  c2->SaveAs("hPolarAngles.root");

  TCanvas *c3 = new TCanvas("c3", "Canvas for hAzimuthalAngles", 800, 600);
  hAzimuthalAngles->Draw();
  c3->SaveAs("hAzimuthalAngles.pdf");
  c3->SaveAs("hAzimuthalAngles.C");
  c3->SaveAs("hAzimuthalAngles.root");

  TCanvas *c4 = new TCanvas("c4", "Canvas for hImpulse", 800, 600);
  hImpulse->Draw();
  c4->SaveAs("hImpulse.pdf");
  c4->SaveAs("hImpulse.C");
  c4->SaveAs("hImpulse.root");

  TCanvas *c5 = new TCanvas("c5", "Canvas for hTransverseImpulse", 800, 600);
  hTransverseImpulse->Draw();
  c5->SaveAs("hTransverseImpulse.pdf");
  c5->SaveAs("hTransverseImpulse.C");
  c5->SaveAs("hTransverseImpulse.root");

  TCanvas *c6 = new TCanvas("c6", "Canvas for hEnergy", 800, 600);
  hEnergy->Draw();

  TLine *line1 = new TLine(0.13957, 0, 0.13957,
                           1E5);  // Disegna la retta verticale a x = 3
  line1->SetLineColor(kRed);      // Imposta il colore della linea
  line1->SetLineWidth(2);         // Imposta lo spessore della linea
  line1->Draw();

  TLine *line2 = new TLine(0.49367, 0, 0.49367,
                           1E5);  // Disegna la retta verticale a x = 3
  line2->SetLineColor(kRed);      // Imposta il colore della linea
  line2->SetLineWidth(2);         // Imposta lo spessore della linea
  line2->Draw();

  TLine *line3 = new TLine(0.93827, 0, 0.93827,
                           1E5);  // Disegna la retta verticale a x = 3
  line3->SetLineColor(kRed);      // Imposta il colore della linea
  line3->SetLineWidth(2);         // Imposta lo spessore della linea
  line3->Draw();

  TLine *line4 = new TLine(0.89166, 0, 0.89166,
                           1E5);  // Disegna la retta verticale a x = 3
  line4->SetLineColor(kRed);      // Imposta il colore della linea
  line4->SetLineWidth(2);         // Imposta lo spessore della linea
  line4->Draw();

  c6->SaveAs("hEnergy.pdf");
  c6->SaveAs("hEnergy.C");
  c6->SaveAs("hEnergy.root");

  TCanvas *c7 = new TCanvas("c7", "Canvas for hInvMassAll", 800, 600);
  hInvMassAll->Draw();
  c7->SaveAs("hInvMassAll.pdf");
  c7->SaveAs("hInvMassAll.C");
  c7->SaveAs("hInvMassAll.root");

  TCanvas *c8 = new TCanvas("c8", "Canvas for hInvMassCon", 800, 600);
  hInvMassCon->Draw();
  c8->SaveAs("hInvMassCon.pdf");
  c8->SaveAs("hInvMassCon.C");
  c8->SaveAs("hInvMassCon.root");

  TCanvas *c9 = new TCanvas("c9", "Canvas for hInvMassDis", 800, 600);
  hInvMassDis->Draw();
  c9->SaveAs("hInvMassDis.pdf");
  c9->SaveAs("hInvMassDis.C");
  c9->SaveAs("hInvMassDis.root");

  TCanvas *c10 = new TCanvas("c10", "Canvas for hInvMassConPiK", 800, 600);
  hInvMassConPiK->Draw();
  c10->SaveAs("hInvMassConPiK.pdf");
  c10->SaveAs("hInvMassConPiK.C");
  c10->SaveAs("hInvMassConPiK.root");

  TCanvas *c11 = new TCanvas("c11", "Canvas for hInvMassDisPiK", 800, 600);
  hInvMassDisPiK->Draw();
  c11->SaveAs("hInvMassDisPiK.pdf");
  c11->SaveAs("hInvMassDisPiK.C");
  c11->SaveAs("hInvMassDisPiK.root");

  TCanvas *c12 = new TCanvas("c12", "Canvas for hBenchmark", 800, 600);
  hBenchmark->Draw();
  c12->SaveAs("hBenchmark.pdf");
  c12->SaveAs("hBenchmark.C");
  c12->SaveAs("hBenchmark.root");

  TCanvas *c13 = new TCanvas("c13", "Canvas for hKStar", 800, 600);
  hKStar->Draw();
  c13->SaveAs("hKStar.pdf");
  c13->SaveAs("hKStar.C");
  c13->SaveAs("hKStar.root");

  TCanvas *c14 = new TCanvas("c14", "Canvas for hDiff1", 800, 600);
  hDiff1->Draw();
  //fitKStar1->Draw("same");
  c14->SaveAs("hDiff1.pdf");
  c14->SaveAs("hDiff1.C");
  c14->SaveAs("hDiff1.root");

  TCanvas *c15 = new TCanvas("c15", "Canvas for hDiff2", 800, 600);
  hDiff2->Draw();
  fitKStar2->Draw("same");
  c15->SaveAs("hDiff2.pdf");
  c15->SaveAs("hDiff2.C");
  c15->SaveAs("hDiff2.root");

  TCanvas *c16 = new TCanvas("c16", "Cosa succede", 800, 600);
  hInvMassDisPiK->Draw();
  hInvMassConPiK->Draw("SAME");
  c16->Update();
  c16->SaveAs("Cosa.pdf");
  c16->SaveAs("Cosa.C");
  c16->SaveAs("Cosa.root");

  TBrowser *browser = new TBrowser();
}
