#include "particle.hpp"

const size_t Particle::MaxNumParticleType_ = 10;
Particle::Particle(std::string NameParticle, double Px = 0, double Py = 0,
                   double Pz = 0)
    : P.x(Px),
P.y(Py), P.z(Pz) {}