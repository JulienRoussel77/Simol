#ifndef SIMOL_SPACESINUS_HPP
#define SIMOL_SPACESINUS_HPP

#include "Potential.hpp"

namespace simol
{
  
  class SpaceSinus : public Potential
  {
  public:
    SpaceSinus(Input const& input);
    double operator()(Vector<double> const& position) const;
    Vector<double> gradient(Vector<double> const& position) const;
    double laplacian(Vector<double> const& position) const;
    
  private:
    double amplitude_;
    double pulsation_;
  };
  
}

#endif