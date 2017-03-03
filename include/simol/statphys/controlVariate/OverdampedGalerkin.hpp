#ifndef SIMOL_OVDGALERKIN_HPP
#define SIMOL_OVDGALERKIN_HPP

#include "Galerkin.hpp"

namespace simol
{
  class OverdampedGalerkin : public Galerkin
  {
  public:
    OverdampedGalerkin(Input const& input);
    virtual void createOperators();
    
    //void createLeta();
    //virtual void computeExpToTrigTens();
    virtual void compute() = 0;
    virtual DVec CVcoeffsVec() const;
    virtual DVec CVObservable() const = 0;
  };
  
  class PeriodicOverdampedGalerkin : public OverdampedGalerkin
  {
  public:
    PeriodicOverdampedGalerkin(Input const& input);
    
    //DVec getGradV() const;
    virtual void compute();
    virtual DVec CVObservable() const;
  };
  
  class ColloidOverdampedGalerkin : public OverdampedGalerkin
  {
  public:
    ColloidOverdampedGalerkin(Input const& input);
    
    virtual void compute();
    virtual DVec CVObservable() const;
  };
}

#endif