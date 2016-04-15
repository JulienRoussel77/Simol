#include "simol/statphys/potential/DoubleWell.hpp"

namespace simol
{
  
  DoubleWell::DoubleWell(Input const & input):Potential(input), height_(input.height()), interWell_(input.interWell())
  {}
  
  double DoubleWell::operator()(double position) const
  { return height_*pow(position-interWell_/2, 2)*pow(position+interWell_/2, 2); }

  Vector<double> DoubleWell::gradient(double position) const
  { 
    return Vector<double>(1, 4*height_*position*(position-interWell_/2)*(position+interWell_/2));
  }
  
}
