#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP

#include <cstdlib>
#include <iostream>

class ParticleType {
 public:
  inline const char* GetName() const { return Name_; }
  inline double GetMass() const { return Mass_; }
  inline int GetCharge() const { return Charge_; }
  virtual void Print() const;
  ParticleType(const char* Name, const double Mass, const int Charge);

  inline virtual double GetWidth() const { return 0; }

 private:
  const char* Name_;
  const double Mass_;
  const int Charge_;
};

#endif
