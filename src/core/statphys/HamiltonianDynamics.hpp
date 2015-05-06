#ifndef SIMOL_HAMILTONIANDYNAMICS_HPP
#define SIMOL_HAMILTONIANDYNAMICS_HPP


namespace simol
{

  template<class ScalarType>
  class HamiltonianDynamics
  {
    public:
      HamiltonianDynamics(ScalarType mass,
                          Potential<ScalarType> const & potential);

      Potential<ScalarType> const & potential() const;

    private:
      ScalarType mass_;
      Potential<ScalarType> potential_;
  };

}

#include "HamiltonianDynamics.ipp"

#endif
