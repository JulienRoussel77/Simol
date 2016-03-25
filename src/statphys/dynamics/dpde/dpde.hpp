#ifndef SIMOL_DPDE_HPP
#define SIMOL_DPDE_HPP

#include "dynamics.hpp" 

namespace simol
{
  
  class DPDE:public UniformStochasticDynamics
  {
    double heatCapacity_;
    double gamma_;
  public:
    DPDE(Input const&  input);
    virtual double& gamma();  // friction de reference \gamma_\star
    virtual double gamma_DPDE(double intEnergy);  // friction dependant de l'energie interne
    virtual double& heatCapacity();
    double sigma() const;
    void energyReinjection(Particle& particle);
  };
  
}

#endif
