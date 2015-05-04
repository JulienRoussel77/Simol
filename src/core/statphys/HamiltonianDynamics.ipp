#ifndef SIMOL_HAMILTONIANDYNAMICS_IPP
#define SIMOL_HAMILTONIANDYNAMICS_IPP

namespace simol
{

  template<class ScalarType> inline
  HamiltonianDynamics<ScalarType>::HamiltonianDynamics(ScalarType mass,
                                                       Potential<ScalarType> const & potential)
  :mass_(mass),potential_(potential)
  {}

  template<class ScalarType> inline
  Potential<ScalarType> const & HamiltonianDynamics<ScalarType>::potential() const
  {
    return potential_;
  }



}

#endif
