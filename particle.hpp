#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <array>
#include "particle_type.hpp"
#include "resonance_type.hpp"

struct Impulse {
  double x;
  double y;
  double z;
};

class Particle {
 public:
  Particle(std::string NameParticle, double Px, double Py, double Pz);

 private:
  static const size_t MaxNumParticleType_;
  static std::array<ParticleType*, MaxNumParticleType_> ParticleType_;
  static int NParticleType_;
  int Index_;
  Impulse P;
  int FindParticle(std::string NameParticle_);
};

#endif