#ifndef SIMOL_HAMILTONDYNAMICS_IPP
#define SIMOL_HAMILTONDYNAMICS_IPP

namespace simol
{

  template<class ScalarType> inline
  HamiltonDynamics<ScalarType>::HamiltonDynamics(ScalarType mass,
                                                 Potential<ScalarType> const & potential)
  :mass_(mass),potential_(potential)
  {}

  template<class ScalarType> inline
  Potential<ScalarType> const & HamiltonDynamics<ScalarType>::potential() const
  { return potential_; }



}

#endif
