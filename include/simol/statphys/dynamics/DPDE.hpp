#ifndef SIMOL_DPDE_HPP
#define SIMOL_DPDE_HPP

#include "LangevinBase.hpp"
#include "simol/statphys/system/NBody.hpp"

namespace simol
{

  class DPDE: public LangevinBase
  {
    double heatCapacity_;
    double kappa_;
    double cutOff_;
    double rejectionCount_;
    double totalCountForRejection_;
    double negativeEnergiesCount_;
    double rejectionRate_;
  public:
    DPDE(Input const&  input);
    virtual string dynamicsName() const {return "DPDE";}
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

    //-- thermal conduction --
    double& kappa();
    const double& kappa() const;
    
    //-- functions pour integration of the fluctuation/dissipation (isolated systems) --
    void energyReinjection(Particle& particle);
    void metropolizedEnergyReinjection(Particle& particle);

     //-- functions for integration of the fluctuation/dissipation (NBody systems) --
    double pairwiseFluctuationDissipation(double v12, double dist, double internalEnergy1, double internalEnergy2, double reduced_mass);
    void secondPartThermalization(Particle& particle);

    //-- auxiliary functions --
    DVec effectiveDrift(Particle& particle);
    void incrementRejection();
    void incrementTotalCountForRejection();
    double chi(double dist);
    void acceptRejectRate(double v12_current, double v12_init, double internalEnergy1, double internalEnergy2, double reduced_mass, double dist);
    
    virtual void getThermo(Output& output) const;
    virtual void specificComputeOutput(Output& output) const;
    
    void fluctuationDissipation(NBody& syst);
    void elementaryFluctuationDissipation(System const& syst, Particle& particle1, Particle& particle2);
    double pairwiseThermalConduction(double dist, double internalEnergy1, double internalEnergy2);
    
    void thermalConduction(NBody& syst);
    void elementaryThermalConduction(System const& syst, Particle& particle1, Particle& particle2);
  };
  
}

#endif
