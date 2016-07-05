#ifndef SIMOL_DPDE_HPP
#define SIMOL_DPDE_HPP

#include "LangevinBase.hpp"

namespace simol
{

  class DPDE: public LangevinBase
  {
    double heatCapacity_;
  public:
    DPDE(Input const&  input);
    virtual void printName() const;
    virtual double gamma_DPDE(double intEnergy);  // friction depending on internal energies
    virtual double& heatCapacity();
    double sigma() const;
 
    //-- functions pour integration of the fluctuation/dissipation --
    void energyReinjection(Particle& particle);
    void metropolizedEnergyReinjection(Particle& particle);

    //-- auxiliary functions --
    Vector<double> effectiveDrift(Particle& particle);
  };
  
}

#endif
