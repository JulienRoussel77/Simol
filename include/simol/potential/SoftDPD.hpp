#ifndef SIMOL_SOFTDPD_HPP
#define SIMOL_SOFTDPD_HPP

#include "Potential.hpp"

namespace simol
{
  
  class SoftDPD : public Potential
  {
  public:
    SoftDPD(Input const& input);
    double operator()(double dist) const;
    double scalarGradient(double dist) const;
  private:
    
    double epsilon_;
    //double cutOffRadius_;
  };

}

#endif
