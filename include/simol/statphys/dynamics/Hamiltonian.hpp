#ifndef SIMOL_HAMILTONIAN_HPP
#define SIMOL_HAMILTONIAN_HPP

#include "Dynamics.hpp"

namespace simol
{
  class Hamiltonian : public Dynamics
  {
    public:
      Hamiltonian(Input const&  input);
      virtual void printName() const;
      virtual void computeThermo(Output& output) const;
      virtual void computeGeneratorOnBasis(CVBasis& cvBasis, System const& syst) const;
  };

}

#endif
