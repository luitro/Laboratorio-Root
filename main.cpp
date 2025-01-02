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

  // Set the seed for the simulation
  gRandom->SetSeed();

  // Declare and define some parameters of the events
  int N = 120;  // Maximum number of particles per event
  Particle
      EventParticles[N];  // Array that stores all the particles in each event
  Double_t phi;           // Azimuthual angle
  Double_t theta;         // Polar angle
  Double_t p;             // Impulse magnitude
  Double_t
      x;  // Random number for generate the particles in different proportions
  int firstempty;  // Track the index of the first empty element in
                   // EventParticles
  double PiK =
      0.13957 +
      0.49367;  // The sum of pion and kaon masses for invartiant mass cts
  Int_t BINS = 1.2 * 1E3;  // Number of bins for invariant mass histograms

  // Define all the histograms
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
  TH1F *h7 =
      new TH1F("h7", "Invariant mass between all particles", BINS, 0, 10);
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
  TH1F *h12 = new TH1F("h12", "Benchmark", 1.5 * 1E3, 0, 7);
  h12->Sumw2();

  // Outer loop: processes a total of 1E5 events. Each event corresponds to a
  // physical process or simulation step.
  for (int i = 0; i < 1E5; ++i) {
    firstempty = 100;  // Reset the "firstempty" variable for each event.
                       // Inner loop: simulates the generation of 100 particles
                       // for the current event.
    for (Int_t j = 0; j < 100; ++j) {
      // Generate random polar and azimuthual angles using uniform distribution
      theta = gRandom->Uniform(0., TMath::Pi());
      phi = gRandom->Uniform(0., 2 * TMath::Pi());

      // Generate random impulse magnitude using exponanial distribution
      p = gRandom->Exp(1.);

      // Set impulse cartesian componets for the particle
      EventParticles[j].SetP(p * sin(theta) * cos(phi),
                             p * sin(theta) * sin(phi), p * cos(theta));

      // Assign particle type: x is extracted from a uniform distribution in the
      // range [0,1]
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
        // Handle the K* decade in a Pi- and a K*
        EventParticles[j].SetIndex("K*");
        EventParticles[firstempty].SetIndex("K+");
        EventParticles[firstempty + 1].SetIndex("Pi-");
        EventParticles[j].Decay2body(EventParticles[firstempty],
                                     EventParticles[firstempty + 1]);

        // Fill the benchmark histograms
        h12->Fill(
            EventParticles[firstempty].InvMass(EventParticles[firstempty + 1]));

        // Update the first empty index
        firstempty += 2;
      } else {
        // Handle the K* decade in a Pi+ and a K-
        EventParticles[j].SetIndex("K*");
        EventParticles[firstempty].SetIndex("K-");
        EventParticles[firstempty + 1].SetIndex("Pi+");
        EventParticles[j].Decay2body(EventParticles[firstempty],
                                     EventParticles[firstempty + 1]);
        h12->Fill(
            EventParticles[firstempty].InvMass(EventParticles[firstempty + 1]));
        firstempty += 2;
      }

      // Fill histograms for single-particle properties
      h1->Fill(EventParticles[j].GetIndex());
      h2->Fill(theta);
      h3->Fill(phi);
      h4->Fill(p);
      h5->Fill(sqrt(pow(p, 2) + pow(EventParticles[j].GetPx(), 2)));
      h6->Fill(EventParticles[j].Energy());
    }

    // Analyze particle pairs for invariant mass histograms
    for (int a = 0; a < firstempty; ++a) {
      // Exclude the first particle to be null or K*
      if (EventParticles[a].GetIndex() != -1 &&
          EventParticles[a].GetIndex() != 6) {
        for (int b = a + 1; b < firstempty; ++b) {
          // Exclude the second particle to be null or K*
          if (EventParticles[b].GetIndex() != -1 &&
              EventParticles[b].GetIndex() != 6) {
            // Calculate the invariant mass
            Double_t c = EventParticles[a].InvMass(EventParticles[b]);
            // Fill the histogram with the invariant mass between all particle
            h7->Fill(c);
            /* Check charge concordance and fill corresponding histograms:
            positive charged particles have an odd index and negative charged
            particles have an even index. The sum of the indexes is even, if
            both are even or odd; instead, the sum of the
            indexes is odd, if one is even and the other one is odd.*/
            if (((EventParticles[a].GetIndex() + EventParticles[b].GetIndex()) %
                 2) == 0) {
              h8->Fill(c);  // Same charge
              if (EventParticles[a].GetMass() + EventParticles[b].GetMass() ==
                  PiK)
                h10->Fill(c);  // Same charge Pi and K
            } else {
              h9->Fill(c);  // Opposite charge
              if (EventParticles[a].GetMass() + EventParticles[b].GetMass() ==
                  PiK)
                h11->Fill(c);  // Opposite charge Pi and K
            }
          }
        }
      }
    }
  }

  // Save histograms to a ROOT file
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
  outputFile->Close();
}