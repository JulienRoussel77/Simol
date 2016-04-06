#ifndef SIMOL_SUMSINUSOIDAL_HPP
#define SIMOL_SUMSINUSOIDAL_HPP

#include "Potential.hpp"

namespace simol
{
  
  class SumSinusoidal : public Potential{
  public:
    SumSinusoidal(Input const& input);
    double operator()(double position) const;
    Vector<double> gradient(double position) const;
    virtual double laplacian(double position) const;
    
    
  private:
    double amplitude_;
    double pulsation_;
  };
  
}

#endif