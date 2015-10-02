#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP


#include <iostream>
#include "particle.hpp"

namespace simol
{

  class Output
  {
    std::string outputFilename_;
  public:
    Output(std::string const& outputFilename);
    void display(double const& time, Particle const& particle);
  };

}
#endif