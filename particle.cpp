#include "particle.hpp"

#include <cmath>

ParticleType* Particle::Particles_[MaxNumParticleType_];
int Particle::NParticleType_ = 0;

Particle::Particle(const char* NameParticle, double Px, double Py, double Pz)
    : Index_(FindParticle(NameParticle)), P_({Px, Py, Pz}) {}

bool Particle::AddParticleType(const char* Name, const double Mass,
                               const int Charge, const double Width) {
  if (FindParticle(Name) != -1) {
    std::cout << "Particle " << Name << " already exists\n";
    return false;
  }
  if (NParticleType_ > static_cast<int>(MaxNumParticleType_)) {
    std::cout << "You have reached the maximum particle type limit";
    return false;
  }
  if (Width == 0) {
    Particles_[NParticleType_] = new ParticleType(Name, Mass, Charge);
    NParticleType_++;
    std::cout << "Particle " << Name << " has been added correctly" << std::endl;
    return true;
  } else {
    Particles_[NParticleType_] = new ResonanceType(Name, Mass, Charge, Width);
    NParticleType_++;
    std::cout << "Particle" << Name << "has been added correctly" << std::endl;
    return true;
  }
  return false;
}

void Particle::SetIndex(const int Index) {
  if (Index != Index_) Index_ = Index;
}

void Particle::SetIndex(const char* Name) {
  if (FindParticle(Name) != -1) Index_ = FindParticle(Name);
}

int Particle::FindParticle(const char* NameParticle_) {
  for (int i = 0; i < NParticleType_; i++) {
    if (Particles_[i] != nullptr && Particles_[i]->GetName() == NameParticle_) {
      return i;
    }
  }
  /* if (NameParticle_ != "Null") {
    std::cout << "Particle " << NameParticle_ << " not found" << std::endl;
  } */ 
  return -1;
}

void Particle::PrintParticles() {
  std::cout << "Particle Table:\n";
  for (const auto& i : Particles_) {
    if (i != nullptr)
      std::cout << i->GetName() << '\t';
    else {
      std::cout << std::endl;
      break;
    }
  }
}

void Particle::PrintParticle() {
  std::cout << "Particle name: " << Particles_[Index_]->GetName()
            << "\nParticle index: " << Index_
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