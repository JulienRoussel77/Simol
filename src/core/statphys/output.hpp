#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP


#include <iostream>
#include "particle.hpp"

namespace simol
{

  class Output
  {
    std::string outputFilename_;
    std::ofstream outAtoms_;
    std::ofstream outReplica_;
  public:
    Output(std::string const& outputFilename);
    void display(Particle const& particle, double time = 0);
    void finalDisplay(Particle const& particle, double time = 0);
  };

}
#endif