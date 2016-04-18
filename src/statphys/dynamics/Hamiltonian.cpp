
#include "simol/statphys/dynamics/Hamiltonian.hpp"

namespace simol
{
  //! Constructs a Hamiltonian dynamics (constant energy)
  Hamiltonian::Hamiltonian(Input const& input)
    : Dynamics(input)
  {}

  void Hamiltonian::printName() const
  {
    std::cout << "DynamicsType = Hamiltonian" << std::endl;
  }
}
