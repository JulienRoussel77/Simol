#ifndef SIMOL_DPDE_HPP
#define SIMOL_DPDE_HPP

#include "LangevinBase.hpp"

namespace simol
{

  class DPDE: public LangevinBase
  {
    double heatCapacity_;
    double cutOff_;
    double rejectionCount_;
    double totalCountForRejection_;
    double rejectionRate_;
    double negativeEnergiesCount_;
  public:
    DPDE(Input const&  input);
    virtual void printName() const;
    virtual double& cutOff();
    double sigma() const;
    const double& rejectionCount() const;
    double& rejectionCount();
    const double& totalCountForRejection() const;
    double& totalCountForRejection();
    const double& rejectionRate() const;
    double& rejectionRate();
    const double& negativeEnergiesCount() const;
    double& negativeEnergiesCount();

    //-- microscopic EOS --
    virtual double gamma_DPDE(double intEnergy);  
    virtual double internalTemperature(double intEnergy) const;
    double& heatCapacity();
    const double& heatCapacity() const;
    
    //-- functions pour integration of the fluctuation/dissipation (isolated systems) --
    void energyReinjection(Particle& particle);
    void metropolizedEnergyReinjection(Particle& particle);

     //-- functions for integration of the fluctuation/dissipation (NBody systems) --
    double pairwiseFluctuationDissipation(double v12, double dist, double internalEnergy1, double internalEnergy2, double reduced_mass);
    void secondPartThermalization(Particle& particle);

    //-- auxiliary functions --
    Vector<double> effectiveDrift(Particle& particle);
    void incrementRejection();
    void incrementTotalCountForRejection();
    double chi(double dist);
    void acceptRejectRate(double v12_current, double v12_init, double internalEnergy1, double internalEnergy2, double reduced_mass, double dist);
    
    virtual void specificComputeOutput(Output& output) const;
  };
  
}

#endif
