#include "particle.hpp"

#include <cmath>

Particle::Particle(const char* NameParticle, double Px = 0, double Py = 0,
                   double Pz = 0)
    : P_({Px, Py, Pz}) {
  Index_ = FindParticle(NameParticle);
}
void Particle::AddParticleType(const char* Name, const double Mass,
                               const int Charge, const double Width) {
  if (NParticleType_ == 0) {
    Particles_[NParticleType_] = &(ParticleType(Name, Mass, Charge));
  }
  if (FindParticle(Name) != -1)
    std::cout << "Particle " << Name << " already exists\n";
  else {
    NParticleType_++;
    Particles_[NParticleType_] = &(ParticleType(Name, Mass, Charge));
  }
};

void Particle::SetIndex(const int Index) {
  if (Index != Index_) Index_ = Index;
}

void Particle::SetIndex(const char* Name) {
  if (FindParticle(Name) != -1) Index_ = FindParticle(Name);
}

int Particle::FindParticle(const char* NameParticle_) {
  for (int i = 0; i < MaxNumParticleType_; i++) {
    if ((*Particles_[i]).GetName() == NameParticle_)
      return i;
    else {
      std::cout << "Particle name not corresponding" << std::endl;
      return -1;
    }
  }
}

void Particle::PrintParticles() {
  for (const auto& i : Particles_) {
    std::cout << "Particle Table:\n" << i->GetName() << '\t';
  }
}

void Particle::PrintParticle() {
  std::cout << "Particle index: " << Index_
            << "\nParticle name: " << Particles_[Index_]->GetName()
            << "\nParticle impulse components: (" << P_.x << ", " << P_.y
            << ", " << P_.z << ")\n";
}

double Particle::Energy() const {
  return sqrt(pow(this->GetMass(), 2) + pow(this->GetPx(), 2) +
              pow(this->GetPy(), 2) + pow(this->GetPz(), 2));
}

double Particle::InvMass(const Particle& p) {
  return sqrt(pow((this->Energy() + p.Energy()), 2) -
              pow(this->GetPx() + p.GetPx(), 2) -
              pow(this->GetPy() + p.GetPy(), 2) -
              pow(this->GetPz() + p.GetPz(), 2));
}

void Particle::SetP(double px, double py, double pz) {
  P_.x = px;
  P_.y = py;
  P_.z = pz; 
}