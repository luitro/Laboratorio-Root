#include "resonance_type.hpp"

void ResonanceType::Print() const {
  ParticleType::Print();
  std::cout << "Width: " << Width_ << std::endl;
}

ResonanceType::ResonanceType(const char* Name, const double Mass, const int Charge, const double Width)
    : ParticleType(Name, Mass, Charge), Width_(Width) {};
