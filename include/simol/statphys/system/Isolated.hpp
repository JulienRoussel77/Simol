#ifndef SIMOL_ISOLATED_HPP
#define SIMOL_ISOLATED_HPP
#include "System.hpp"

namespace simol
{

  class Isolated : public System
  {
    public:
      Isolated(Input const& input);
      void printName() const;
      void computeAllForces();
      virtual void thermalize(Dynamics& model);
      virtual void computeKineticEnergy(Output& output) const;
      virtual void computePotentialEnergy(Output& output) const;
      //virtual void computePressure(Output& output) const;
      //virtual void computeInternalEnergy(Output& output) const;
      //virtual void computeInternalTemperature(Output& output) const;
  };


}

#endif