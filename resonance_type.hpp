#ifndef RESONANCE_TYPE_HPP
#define RESONANCE_TYPE_HPP

#include "particle_type.hpp"

class ResonanceType : public ParticleType {
 public:
  inline const double GetWidth() { return Width_; }
  void Print() const override;
  ResonanceType(char Name, double Mass, int Charge, double Width);

 private:
  const double Width_;
};

#endif