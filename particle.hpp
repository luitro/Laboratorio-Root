#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "resonance_type.hpp"

struct Impulse {
  double x;
  double y;
  double z;
};

class Particle {
 public:
  Particle(const char* NameParticle, double Px = 0, double Py = 0, double Pz = 0);
  inline int GetIndex() const { return Index_; }
  static bool AddParticleType(const char* Name, const double Mass,
                              const int Charge, const double Width = 0);
  void SetIndex(const int Index);
  void SetIndex(const char* Name);
  static void PrintParticles();
  void PrintParticle();
  inline double GetPx() const { return P_.x; }
  inline double GetPy() const { return P_.y; }
  inline double GetPz() const { return P_.z; }
  inline double GetMass() const { return Particles_[Index_]->GetMass(); }
  double Energy() const; 
  double InvMass(const Particle& p);
  void SetP(double px, double py, double pz);

  int Decay2body(Particle &dau1,Particle &dau2) const;

 private:
  static const size_t MaxNumParticleType_ = 10;
  static ParticleType* Particles_[MaxNumParticleType_];
  static int NParticleType_;
  int Index_;
  Impulse P_;
  static int FindParticle(const char* NameParticle_);

  void Boost(double bx, double by, double bz);
};

#endif