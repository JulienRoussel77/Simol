#ifndef SIMOL_POTENTIAL_HPP
#define SIMOL_POTENTIAL_HPP

namespace simol
{

  template<class ScalarType>
  class Potential
  {
    public:
      Potential(ScalarType parameter, ScalarType pulsatance);
      ScalarType operator()(ScalarType const & position) const;
      ScalarType derivative(ScalarType const & position) const;
      ScalarType force(ScalarType const & position) const;

    private:
      ScalarType parameter_;
      ScalarType energy_;
  };

}

#include "Potential.ipp"

#endif
