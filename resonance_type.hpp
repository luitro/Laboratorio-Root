#ifndef RESONANCE_TYPE_HPP
#define RESONANCE_TYPE_HPP

#include "particle_type.hpp"

class ResonanceType : public ParticleType {
 public:
  inline double GetWidth() const override { return Width_; }
  void Print() const override;
  ResonanceType(const char* Name, const double Mass, const int Charge,
                const double Width);

 private:
  const double Width_;
};

#endif