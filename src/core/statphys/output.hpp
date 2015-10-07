#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP


#include <iostream>
#include "particle.hpp"

namespace simol
{

  class Output
  {
    std::string outputFilename_;
    std::ofstream out_;
  public:
    Output(std::string const& outputFilename);
    void display(Particle const& particle, double time = 0);
  };

}
#endif