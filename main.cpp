#include "particle.hpp"

int main() {
  Particle::AddParticleType("P+", 1., 1 );

  Particle Protone("P+");
  Protone.PrintParticle(); 

  Particle::AddParticleType("P-", 1., -1);
  Particle coso("P-");
  coso.PrintParticle(); 
  Particle::PrintParticles();
}