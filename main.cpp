#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"
#include "particle.hpp"

int main() {
  Particle::AddParticleType("Pi+", 0.13957, 1);
  Particle::AddParticleType("Pi-", 0.13957, -1);
  Particle::AddParticleType("K+", 0.49367, 1);
  Particle::AddParticleType("K-", 0.49367, -1);
  Particle::AddParticleType("P+", 0.93827, 1);
  Particle::AddParticleType("P-", 0.93827, -1);
  Particle::AddParticleType("K*", 0.89166, 0, 0.050);

  gRandom->SetSeed();

  int N = 120;
  Particle EventParticles[N] = 0;
  Double_t phi;
  Double_t theta;
  Double_t p;
  Double_t x;
  int firstempty = 100;

  TH1F *h1 = new TH1F("h1", "Particle types", 7, 0, 1E7);
  TH1F *h2 = new TH1F("h2", "Polar angles", 1E3, 0, 2 * TMath::Pi());
  TH1F *h3 = new TH1F("h3", "Azimuthal angles", 1E3, 0, TMath::Pi());
  TH1F *h4 = new TH1F("h4", "Impulse", 1E4, 0, 10);
  TH1F *h5 = new TH1F("h5", "Transverse impulse", 1E4, 0, 10);
  TH1F *h6 = new TH1F("h6", "Energy", 1E4, 0, 10);
  TH1F *h7 = new TH1F("h7", "Invariant mass", 1E5, 0, 10);

  for (int i = 0; i < 1E5; ++i) {
    for (Int_t j = 0; j < 100; ++i) {
      phi = gRandom->Uniform(0., 2 * TMath::Pi());
      theta = gRandom->Uniform(0., TMath::Pi());
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
        firstempty += 2;
      } else {
        EventParticles[j].SetIndex("K*");
        EventParticles[firstempty].SetIndex("K-");
        EventParticles[firstempty + 1].SetIndex("Pi+");
        EventParticles[j].Decay2body(EventParticles[firstempty],
                                     EventParticles[firstempty + 1]);
        firstempty += 2;
      }

      h1->Fill(EventParticles[j].GetIndex());
      h2->Fill(theta);
      h3->Fill(phi);
      h4->Fill(p);
      h5->Fill(sqrt(pow(p, 2) + pow(EventParticles[j].GetPx(), 2)));
      h6->Fill(EventParticles[j].Energy());
    }
    for (auto& const i : EventParticles){
      h7->Fill(i.InvMass())
    }
  }
}