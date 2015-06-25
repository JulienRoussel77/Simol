#ifndef SIMOL_HAMILTONDYNAMICS_HPP
#define SIMOL_HAMILTONDYNAMICS_HPP


namespace simol
{

  template<class ScalarType>
  class HamiltonDynamics
  {
    public:
      HamiltonDynamics(Potential<ScalarType> const & potential);

      Potential<ScalarType> const & potential() const;

    private:
      //ScalarType mass_;
      Potential<ScalarType> potential_;
  };

}

#include "HamiltonDynamics.ipp"

#endif
