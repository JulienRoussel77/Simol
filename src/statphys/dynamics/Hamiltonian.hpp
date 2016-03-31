#ifndef SIMOL_HAMILTONIAN_HPP
#define SIMOL_HAMILTONIAN_HPP

#include "Dynamics.hpp"

namespace simol
{
  class Hamiltonian : public Dynamics
  {
    public:
        Hamiltonian(Input const&  input);
  };

}

#endif
