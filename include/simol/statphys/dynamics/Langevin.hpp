#ifndef SIMOL_LANGEVIN_HPP
#define SIMOL_LANGEVIN_HPP

#include "LangevinBase.hpp"

namespace simol
{
  class Langevin : public LangevinBase
  {
    public:
      Langevin(Input const& input);
      virtual string dynamicsName() const {return "Langevin";}
      
      double sigma() const;
      //virtual void updateAfter(Particle& particle);
      virtual void computeGeneratorOnBasis(shared_ptr<CVBasis> cvBasis, System const& syst) const;
    protected:
      double sigma_;
  };
  
  
  class ConstrainedLangevin : public Langevin
  {
    public:
      ConstrainedLangevin(Input const& input);
      virtual double& lagrangeMultiplier() {return lagrangeMultiplier_;}
      virtual const double& lagrangeMultiplier() const {return lagrangeMultiplier_;}
      virtual double& drift() {return drift_;}
      virtual const double& drift() const {return drift_;}
      virtual string dynamicsName() const {return "ConstrainedLangevin";}
    protected:
      double lagrangeMultiplier_;
      double drift_;
  };

}

#endif
