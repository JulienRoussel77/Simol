#ifndef POTENTIAL_IMPL_HPP
#define POTENTIAL_IMPL_HPP

#include "potential.hpp"

#include <cmath>

using std::cout; 
using std::endl; 

namespace simol
{

  Potential::Potential(){}
  
  Potential::Potential(Input const & input, int const& indexOfReplica){}
  
  Potential* createPotential(Input const& input, int const& indexOfReplica)
  {
    if (input.potentialName() == "Sinusoidal")
      return new Sinusoidal(input);
    else if (input.potentialName() == "SumSinusoidal")
      return new SumSinusoidal(input);
    else if (input.potentialName() == "FracSinusoidal")
      return new FracSinusoidal(input);
    else if (input.potentialName() == "DoubleWell")
      return new DoubleWell(input);
    else if (input.potentialName() == "HarmonicWell")
      return new HarmonicWell(input);
    else if (input.potentialName() == "Harmonic")
      return new Harmonic(input);
    else
      std::cout << input.potentialName() << " is not a valid potential !" << std::endl;
    
    return 0;
  }
  
  dvec Potential::force(dvec const & position) const
  { return -derivative(position); }
  
//#### Sinusoidal #####
  
  Sinusoidal::Sinusoidal(Input const & input, int const& indexOfReplica):Potential(input, indexOfReplica), amplitude_(input.amplitude()), pulsation_(2*M_PI/input.length())
  {}
  
  double Sinusoidal::operator()(dvec const & position) const
  { return amplitude_* (1-std::cos(pulsation_*position(0))); }

  dvec Sinusoidal::derivative(dvec const & position) const
  { 
    dvec deriv(1);
    deriv(0) = amplitude_*pulsation_*std::sin(pulsation_*position(0));
    return deriv;
  }
  
  double Sinusoidal::laplacian(dvec const & position) const
  {
    double q = position(0);
    return -amplitude_ * pow(pulsation_, 2) * sin(pulsation_ * q);
  }
  
    
//#### SumSinusoidal #####
  
  SumSinusoidal::SumSinusoidal(Input const & input, int const& indexOfReplica):
    Potential(input, indexOfReplica), 
    amplitude_(input.amplitude()), 
    pulsation_(2*M_PI/input.length())
  {}
  
  double SumSinusoidal::operator()(dvec const & position) const  { return amplitude_* (std::sin(pulsation_*position(0)) 
    
	+ cos(2 * pulsation_*position(0)) 
	+ cos(3 * pulsation_*position(0)) / 3); }

  dvec SumSinusoidal::derivative(dvec const & position) const
  { 
    dvec deriv(1);
    deriv(0) = amplitude_*pulsation_*(cos(pulsation_*position(0))
	  - 2 * sin(2 * pulsation_*position(0))
	  - sin(3 * pulsation_*position(0)));
    return deriv;
  }

    
  double SumSinusoidal::laplacian(dvec const & position) const
  {
    double q = position(0);
    return amplitude_ * pow(pulsation_, 2) * (-sin(pulsation_ * q)
	      - 4 * cos(2 * pulsation_ * q)
	      - 3 * cos(3 * pulsation_ * q));
  }
  
  
  //#### FracSinusoidal #####
  
  FracSinusoidal::FracSinusoidal(Input const & input, int const& indexOfReplica):
    Potential(input, indexOfReplica),
    amplitude_(input.amplitude()),
    pulsation_(2*M_PI/input.length())
  {}
  
  double FracSinusoidal::operator()(dvec const & position) const
  { 
    double q = position(0);
    return amplitude_ * cos(2 * M_PI * q) / (2 + sin(2 * M_PI * q));
  }

  dvec FracSinusoidal::derivative(dvec const & position) const
  { 
    double q = position(0);
    dvec deriv(1);
    //deriv(0) = -(2 * M_PI * sin(2 * M_PI * q))/(2 + sin(2 * M_PI * q))
	//	-(2 * M_PI * pow(cos(2 M_PI * q),2))/(sin(2 * M_PI * q)+2)^2;
    deriv(0) = -amplitude_ * (4 * M_PI * sin(2 * M_PI * q) + 2 * M_PI)/pow(sin(2 * M_PI * q) + 2, 2);
    return deriv;
  }

    
  double FracSinusoidal::laplacian(dvec const & position) const
  {
    double q = position(0);    
    return -amplitude_ * 32 * pow(M_PI,2) * pow(sin(M_PI/4-M_PI * q), 3) * sin(M_PI * q + M_PI/4)
	  / pow(sin(2 * M_PI * q) + 2, 3);
  }
  
//#### DoubleWell #####  
  
    DoubleWell::DoubleWell(Input const & input, int const& indexOfReplica):Potential(input, indexOfReplica), height_(input.height()), interWell_(input.interWell())
  {}
  
  double DoubleWell::operator()(dvec const & position) const
  { return height_*pow(position(0)-interWell_/2, 2)*pow(position(0)+interWell_/2, 2); }

  dvec DoubleWell::derivative(dvec const & position) const
  { 
    dvec deriv(1);
    deriv(0) = 4*height_*position(0)*(position(0)-interWell_/2)*(position(0)+interWell_/2);
    return deriv;
  }
  
  //#### HarmonicWell #####  
  
  HarmonicWell::HarmonicWell(Input const & input, int const& indexOfReplica):Potential(input, indexOfReplica), stiffness_(input.stiffness())
  {}
  
  
  double HarmonicWell::operator()(dvec const & position) const
  { return stiffness_ * pow(position(0), 2); }

  dvec HarmonicWell::derivative(dvec const & position) const
  { 
    dvec deriv(1);
    deriv(0) = 2 * stiffness_ * position(0);
    return deriv;
  }

  
//#### Harmonic #####  
  
  Harmonic::Harmonic(Input const & input, int const& indexOfReplica):Potential(input, indexOfReplica), stiffness_(input.stiffness())
  {}
  
  
  double Harmonic::operator()(dvec const & vecDistance) const
  { return stiffness_ * pow(vecDistance(0) - 1, 2); }

  dvec Harmonic::derivative(dvec const & vecDistance) const
  { 
    dvec deriv(1);
    deriv(0) = 2 * stiffness_ * (1 - vecDistance(0));
    return deriv;
  }

  
}


#endif
