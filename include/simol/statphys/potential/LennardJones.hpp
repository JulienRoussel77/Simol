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
      DVec gradient(double dist) const;
    protected:
      double splineFunction(double reducedDist) const;
      double splineFunctionDerivative(double dist) const;
      double untruncated(double dist) const;
      double untruncatedDerivative(double dist) const;

      double epsilon_;
      double sigma_;
      double cutOffRadius_;
      double splineRatio_;
      double splineRadius_;
      double A_spline_;
      double B_spline_;
      double C3_spline_;
      double C4_spline_;
  };
  
  class WCA : public LennardJones
  {
  public:
    WCA(Input const& input);
    double operator()(double dist) const;
    DVec gradient(double dist) const;
  };

}

#endif