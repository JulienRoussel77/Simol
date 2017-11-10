#ifndef SIMOL_COULOMB_HPP
#define SIMOL_COULOMB_HPP

#include "Potential.hpp"

namespace simol
{
  
  class Coulomb : public Potential
  {
  public:
    Coulomb(Input const& input);
    double operator()(double dist) const;
    double scalarGradient(double dist) const;
  private:
    
    double epsilon_;
    double sigma_;
    double cutOffRadius_;
    double coeff_;
  };

}

#endif
