#include "particle.hpp"

Particle::Particle(const char* NameParticle, double Px = 0, double Py = 0,
                   double Pz = 0)
    : P_({Px,Py,Pz}) {
        Index_ = FindParticle(NameParticle);
    }

int Particle::FindParticle(const char* NameParticle_) {
    for (int i = 0; i < MaxNumParticleType_; i++) {
      if ((*Particles_[i]).GetName() == NameParticle_)
        return i;
      else std::cout << "Particle name not corresponding" << std::endl; 
    }
  }