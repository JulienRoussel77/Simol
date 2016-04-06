#ifndef SIMOL_LENNARDJONES_HPP
#define SIMOL_LENNARDJONES_HPP

#include "Potential.hpp"

namespace simol
{
  
  class LennardJones : public Potential
  {
  public:
    LennardJones(Input const& input);
    double operator()(double dist) const;
    Vector<double> gradient(double dist) const;
  private:
    double epsilon_;
    double sigma_;
    double cutOffRadius_;
    double splineRadius_;
    double splineFunction(double reducedDist) const;
    double untruncated(double dist) const;
    double A_spline_;
    double B_spline_;
  };
  
}

#endif