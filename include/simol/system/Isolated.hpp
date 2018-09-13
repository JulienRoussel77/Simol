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
      virtual string name() const {return "Isolated";}

      virtual void samplePositions(DynamicsParameters const& dynaPara);
      void computeAllForces();
  };


}

#endif
