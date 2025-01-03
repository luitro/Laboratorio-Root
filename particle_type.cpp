#include "particle_type.hpp"

void ParticleType::Print() const {
  std::cout << "Name: " << *Name_ << " ;\nMass: " << Mass_
            << " ;\nCharge: " << Charge_ << std::endl;
}

ParticleType::ParticleType(const char* Name, const double Mass, const int Charge)
    : Name_(Name), Mass_(Mass), Charge_(Charge) {}

