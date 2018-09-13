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
      virtual void simulate (System& syst) const;
      
      double sigma() const;
      virtual const double& drift() const {return parameters_.drift();}
    protected:
  };
  
  
  class ConstrainedLangevin : public Langevin
  {
    public:
      ConstrainedLangevin(Input const& input);
      virtual void simulate (System& syst) const;      
    protected:
  };

}

#endif
