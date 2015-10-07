#ifndef POTENTIAL_IMPL_HPP
#define POTENTIAL_IMPL_HPP

#include "potential.hpp"

#include <cmath>

using std::cout; 
using std::endl; 

namespace simol
{

  Potential::Potential(){}
  
  Potential::Potential(Input const & input){}
  
  Potential* createPotential(Input const& input)
  {
    if (input.potentialName() == "Sinusoidal")
      return new Sinusoidal(input);
    else if (input.potentialName() == "DoubleWell")
      return new DoubleWell(input);
    else if (input.potentialName() == "Harmonic")
      return new Harmonic(input);
    else
      std::cout << input.potentialName() << " is not a valid potential !" << std::endl;
    
    return 0;
  }
  
//#### Sinusoidal #####
  
  Sinusoidal::Sinusoidal(Input const & input):amplitude_(input.amplitude()), pulsation_(2*M_PI/input.length())
  {}
  
  double Sinusoidal::operator()(dvec const & position) const
  { return amplitude_* (1-std::cos(pulsation_*position(0))); }

  dvec Sinusoidal::derivative(dvec const & position) const
  { 
    dvec deriv(1);
    deriv(0) = amplitude_*pulsation_*std::sin(pulsation_*position(0));
    return deriv;
  }

  dvec Sinusoidal::force(dvec const & position) const
  { return -derivative(position); }
  
//#### DoubleWell #####  
  
    DoubleWell::DoubleWell(Input const & input):height_(input.height()), interWell_(input.interWell())
  {}
  
  double DoubleWell::operator()(dvec const & position) const
  { return height_*pow(position(0)-interWell_/2, 2)*pow(position(0)+interWell_/2, 2); }

  dvec DoubleWell::derivative(dvec const & position) const
  { 
    dvec deriv(1);
    deriv(0) = 4*height_*position(0)*(position(0)-interWell_/2)*(position(0)+interWell_/2);
    return deriv;
  }

  dvec DoubleWell::force(dvec const & position) const
  { return -derivative(position); }
  
//#### Harmonic #####  
  
  Harmonic::Harmonic(Input const & input):stiffness_(input.stiffness())
  {}
  
  
  double Harmonic::operator()(double const & distance) const
  { return pow(distance - 1, 2); }

  dvec Harmonic::derivative(dvec const & vecDistance) const
  { 
    dvec deriv(1);
    deriv(0) = 2 * (vecDistance(0) - 1);
    return deriv;
  }

  dvec Harmonic::force(dvec const & position) const
  { return -derivative(position); }
  

  
}


#endif
