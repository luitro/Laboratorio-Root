#include "resonance_type.hpp"

void ResonanceType::Print() const {
  ParticleType::Print();
  std::cout << "Width: " << Width_ << std::endl;
}

ResonanceType::ResonanceType(char Name, double Mass, int Charge, double Width)
    : ParticleType(Name, Mass, Charge), Width_(Width) {};
