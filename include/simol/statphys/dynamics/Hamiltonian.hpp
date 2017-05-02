#ifndef SIMOL_HAMILTONIAN_HPP
#define SIMOL_HAMILTONIAN_HPP

#include "Dynamics.hpp"

namespace simol
{
  class Hamiltonian : public Dynamics
  {
    public:
      Hamiltonian(Input const&  input);
      virtual string dynamicsName() const {return "Hamiltonian";}
      virtual void computeThermo(Output& output) const;
      virtual void computeGeneratorOnBasis(shared_ptr<CVBasis> cvBasis, System const& syst) const;
  };

}

#endif
