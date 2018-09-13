#ifndef SIMOL_BOUNDARYLANGEVIN_HPP
#define SIMOL_BOUNDARYLANGEVIN_HPP

#include "LangevinBase.hpp"

namespace simol
{
  class BoundaryLangevin : public LangevinBase
  {
    public:
      BoundaryLangevin(Input const& input);
      string dynamicsName() const {return "BoundaryLangevin";}
      
      virtual const double& betaLeft() const;
      virtual const double& betaRight() const;
      virtual const double& temperatureLeft() const;
      virtual const double& temperatureRight() const;
      double sigmaLeft() const;
      double sigmaRight() const;
      const double& tauBending() const;
      
      virtual void thermalize(System& syst) const;
      virtual void simulate (System& syst) const;

      void initializeMomenta(vector<Particle>& configuration) const;
      virtual void bending(Particle& particle1, Particle& particle2) const;
      
      virtual void computeProfileChain(Output& output, System const& syst, long int iOfStep) const;
    protected:
      /*double deltaTemperature_;
      double temperatureLeft_;
      double temperatureRight_;
      double betaLeft_;
      double betaRight_;
      double tauBending_;*/
  };
  
  
  class ConstrainedBoundaryLangevin : public BoundaryLangevin
  {
    public:
      ConstrainedBoundaryLangevin(Input const& input);
      virtual void simulate (System& syst) const;
      //virtual double& lagrangeMultiplier() const {return lagrangeMultiplier_;}
      //virtual const double& lagrangeMultiplier() const {return lagrangeMultiplier_;}
      virtual const double& flux() const {return parameters_.flux();}
      //virtual string dynamicsName() const {return "ConstrainedBoundaryLangevin";}
    protected:
      
      //double flux_;
  };

}

#endif
