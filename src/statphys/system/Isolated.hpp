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
    void computeAllForces(Dynamics const& model);
    //void computeFinalOutput(Output& /*output*/, Dynamics const& /*dyna*/);
    //void writeFinalOutput(Output& output, Dynamics const& model);
  };


}

#endif