#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP

#include <cstdlib>
#include <iostream>

class ParticleType {
 public:
  inline const char* GetName() const { return Name_; }
  inline const double GetMass() const { return Mass_; }
  inline const int GetCharge() const { return Charge_; }
  virtual void Print() const;
  ParticleType(const char* Name, const double Mass, const int Charge);

 private:
  const char* Name_;
  const double Mass_;
  const int Charge_;
};

#endif
