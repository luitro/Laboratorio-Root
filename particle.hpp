#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "particle_type.hpp"
#include "resonance_type.hpp"

struct Impulse {
  double x;
  double y;
  double z;
};

class Particle {
 public:
  Particle(const char* NameParticle, double Px, double Py, double Pz);

 private:
  static const size_t MaxNumParticleType_ = 10;
  static ParticleType* Particles_[MaxNumParticleType_];
  static int NParticleType_;
  int Index_;
  Impulse P_;
  int FindParticle(const char* NameParticle_);
};

#endif