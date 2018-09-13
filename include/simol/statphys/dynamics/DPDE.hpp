#ifndef SIMOL_DPDE_HPP
#define SIMOL_DPDE_HPP

#include "LangevinBase.hpp"
#include "simol/statphys/system/NBody.hpp"

namespace simol
{

  class DPDE: public LangevinBase
  {
    double heatCapacity_;
    double heatCapacityEinstein_;
    double einsteinTemperature_;
    double kappa_;
    double cutOff_;
    double rejectionRate_;
    bool doMetropolis_;
    bool doProjectionDPDE_;
    int MTSfrequency_;
    double initialInternalTemperature_; 
    //-- counting rejections for fluctuation/dissipation part --
    double rejectionCountFD_;
    double totalCountForRejectionFD_;
    double negativeEnergiesCountFD_;
    //-- counting rejections for thermal conduction --
    double rejectionCountThermal_;
    double totalCountForRejectionThermal_;
    double negativeEnergiesCountThermal_;
  public:
    DPDE(Input const&  input);
    virtual string dynamicsName() const {return "DPDE";}
    virtual double& cutOff();
    double sigma() const;
    const double& rejectionRate() const;
    double& rejectionRate();
    const double& rejectionCountFD() const;
    double& rejectionCountFD();
    const double& totalCountForRejectionFD() const;
    double& totalCountForRejectionFD();
    const double& negativeEnergiesCountFD() const;
    double& negativeEnergiesCountFD();
    const double& rejectionCountThermal() const;
    double& rejectionCountThermal();
    const double& totalCountForRejectionThermal() const;
    double& totalCountForRejectionThermal();
    const double& negativeEnergiesCountThermal() const;
    double& negativeEnergiesCountThermal();
    const int& MTSfrequency() const;
    int& MTSfrequency();
    const bool& doProjectionDPDE() const;
    bool& doProjectionDPDE();
    
    virtual void thermalize(System& syst);
    virtual void simulate (System& syst);

    //-- microscopic EOS --
    virtual double gamma_DPDE(double intEnergy);  
    virtual double internalTemperature(double intEnergy) const;
    virtual double entropy(double intEnergy) const; 
    virtual double entropy_derivative(double intEnergy) const; 
    double& heatCapacity();
    const double& heatCapacity() const;
    double& heatCapacityEinstein();
    const double& heatCapacityEinstein() const;
    double& einsteinTemperature();
    const double& einsteinTemperature() const;
    
    //-- thermal conduction --
    double& kappa();
    const double& kappa() const;
    
    //-- functions pour integration of the fluctuation/dissipation (isolated systems) --
    void energyReinjection(Particle& particle);
    void metropolizedEnergyReinjection(Particle& particle);

     //-- functions for integration of the fluctuation/dissipation (NBody systems) --
    double pairwiseFluctuationDissipation(double v12, double dist, double internalEnergy1, double internalEnergy2, double reduced_mass);
    void Thermalization(Particle& particle);

    //-- auxiliary functions --
    DVec effectiveDrift(Particle& particle);
    void incrementRejectionFD();
    void incrementTotalCountForRejectionFD();
    void incrementRejectionThermal();
    void incrementTotalCountForRejectionThermal();
    double chi(double dist);
    void acceptRejectRate(double v12_current, double v12_init, double internalEnergy1, double internalEnergy2, double reduced_mass, double dist);
    
    virtual void getThermo(Output& output) const;
    virtual void specificComputeOutput(Output& output) const;
    
    void fluctuationDissipation(System& syst);
    void elementaryFluctuationDissipation(System const& syst, Particle& particle1, Particle& particle2);
    double pairwiseThermalConduction(double chi_n, double internalEnergy1, double internalEnergy2);
    
    void thermalConduction(System& syst);
    void elementaryThermalConduction(System const& syst, Particle& particle1, Particle& particle2);
  };
  
}

#endif
