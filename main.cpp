#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom.h>

#include "particle.hpp"

int main() {
  // Create the static table with all the particles type
  Particle::AddParticleType("Pi+", 0.13957, 1, 0);
  Particle::AddParticleType("Pi-", 0.13957, -1, 0);
  Particle::AddParticleType("K+", 0.49367, 1, 0);
  Particle::AddParticleType("K-", 0.49367, -1, 0);
  Particle::AddParticleType("P+", 0.93827, 1, 0);
  Particle::AddParticleType("P-", 0.93827, -1, 0);
  Particle::AddParticleType("K*", 0.89166, 0, 0.050);

  // Set the seed for the 
  gRandom->SetSeed();

  int N = 120;
  Particle EventParticles[N];
  Double_t phi;
  Double_t theta;
  Double_t p;
  Double_t x;
  int firstempty;
  double PiK = 0.13957 + 0.49367;

  Int_t BINS = 1E3;

  TH1F *h1 = new TH1F("h1", "Particle types", 7, 0, 7);
  h1->Sumw2();
  TH1F *h2 = new TH1F("h2", "Polar angles", 1E3, 0, TMath::Pi());
  h2->Sumw2();
  TH1F *h3 = new TH1F("h3", "Azimuthal angles", 1E3, 0, 2 * TMath::Pi());
  h3->Sumw2();
  TH1F *h4 = new TH1F("h4", "Impulse", 1E4, 0, 10);
  h4->Sumw2();
  TH1F *h5 = new TH1F("h5", "Transverse impulse", 1E4, 0, 10);
  h5->Sumw2();
  TH1F *h6 = new TH1F("h6", "Energy", 1E4, 0, 10);
  h6->Sumw2();
  TH1F *h7 = new TH1F("h7", "Invariant mass between all particles", 1E5, 0, 10);
  h7->Sumw2();
  TH1F *h8 = new TH1F(
      "h8", "Invariant mass between concordant charge particles", BINS, 0, 7);
  h8->Sumw2();
  TH1F *h9 =
      new TH1F("h9", "Invariant mass between disconcordant charge particles",
               BINS, 0, 7);
  h9->Sumw2();
  TH1F *h10 = new TH1F(
      "h10", "Invariant mass between concordant charge Pi and K particles",
      BINS, 0, 7);
  h10->Sumw2();
  TH1F *h11 = new TH1F(
      "h11", "Invariant mass between disconcordant charge Pi and K particles",
      BINS, 0, 7);
  h11->Sumw2();
  TH1F *h12 = new TH1F("h12", "Benchmark", BINS, 0, 7);
  h12->Sumw2();
  TH1F *h13 = new TH1F("h13", "K*", BINS, 0, 7);
  h13->Sumw2();

  for (int i = 0; i < 1E5; ++i) {
    firstempty = 100;
    for (Int_t j = 0; j < 100; ++j) {
      theta = gRandom->Uniform(0., TMath::Pi());
      phi = gRandom->Uniform(0., 2 * TMath::Pi());
      p = gRandom->Exp(1.);
      EventParticles[j].SetP(p * sin(theta) * cos(phi),
                             p * sin(theta) * sin(phi), p * cos(theta));
      x = gRandom->Rndm();
      if (x < 0.4)
        EventParticles[j].SetIndex("Pi+");
      else if (x < 0.8)
        EventParticles[j].SetIndex("Pi-");
      else if (x < 0.85)
        EventParticles[j].SetIndex("K+");
      else if (x < 0.90)
        EventParticles[j].SetIndex("K-");
      else if (x < 0.945)
        EventParticles[j].SetIndex("P+");
      else if (x < 0.99)
        EventParticles[j].SetIndex("P-");
      else if (x < 0.995) {
        EventParticles[j].SetIndex("K*");
        EventParticles[firstempty].SetIndex("K+");
        EventParticles[firstempty + 1].SetIndex("Pi-");
        EventParticles[j].Decay2body(EventParticles[firstempty],
                                     EventParticles[firstempty + 1]);
        h12->Fill(
            EventParticles[firstempty].InvMass(EventParticles[firstempty + 1]));
        firstempty += 2;
      } else {
        EventParticles[j].SetIndex("K*");
        EventParticles[firstempty].SetIndex("K-");
        EventParticles[firstempty + 1].SetIndex("Pi+");
        EventParticles[j].Decay2body(EventParticles[firstempty],
                                     EventParticles[firstempty + 1]);
        h12->Fill(
            EventParticles[firstempty].InvMass(EventParticles[firstempty + 1]));
        firstempty += 2;
      }

      h1->Fill(EventParticles[j].GetIndex());
      h2->Fill(theta);
      h3->Fill(phi);
      h4->Fill(p);
      h5->Fill(sqrt(pow(p, 2) + pow(EventParticles[j].GetPx(), 2)));
      h6->Fill(EventParticles[j].Energy());
    }

    for (int a = 0; a < firstempty; ++a) {
      if (EventParticles[a].GetIndex() != -1 &&
          EventParticles[a].GetIndex() != 6) {
        for (int b = a + 1; b < firstempty; ++b) {
          if (EventParticles[b].GetIndex() != -1 &&
              EventParticles[b].GetIndex() != 6) {
            Double_t c = EventParticles[a].InvMass(EventParticles[b]);
            h7->Fill(c);
            if (((EventParticles[a].GetIndex() + EventParticles[b].GetIndex()) %
                 2) == 0) {
              h8->Fill(c);
              if (EventParticles[a].GetMass() + EventParticles[b].GetMass() ==
                  PiK)
                h10->Fill(c);
            } else {
              h9->Fill(c);
              if (EventParticles[a].GetMass() + EventParticles[b].GetMass() ==
                  PiK)
                h11->Fill(c);
            }
          }
        }
      }
    }
  }
  h13->Add(h9, h8, 1, -1);

  TFile *outputFile = new TFile("histograms.root", "RECREATE");
  h1->Write();
  h2->Write();
  h3->Write();
  h4->Write();
  h5->Write();
  h6->Write();
  h7->Write();
  h8->Write();
  h9->Write();
  h10->Write();
  h11->Write();
  h12->Write();
  h13->Write();
  outputFile->Close();
}