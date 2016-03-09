#ifndef POTENTIAL_IMPL_HPP
#define POTENTIAL_IMPL_HPP

#include "potential.hpp"

#include <cmath>

using std::cout; 
using std::endl; 

namespace simol
{

  Potential::Potential(){}
  
  Potential::Potential(Input const & /*input*/, int const& /*indexOfReplica*/){}
  
  Potential* createPotential(Input const& input, int const& /*indexOfReplica*/)
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
    else if (input.potentialName() == "Rotor")
      return new Rotor(input);
    else if (input.potentialName() == "Quadratic")
      return new Quadratic(input);
    else
      std::cout << input.potentialName() << " is not a valid potential !" << std::endl;
    
    return 0;
  }
  
  Vector<double> Potential::force(Vector<double> const & position) const
  { return -derivative(position); }
  
//#### Sinusoidal #####
  
  Sinusoidal::Sinusoidal(Input const & input, int const& indexOfReplica):Potential(input, indexOfReplica), amplitude_(input.amplitude()), pulsation_(2*M_PI/input.length())
  {}
  
  double Sinusoidal::operator()(Vector<double> const & position) const
  { 
		return amplitude_* (1-std::cos(pulsation_*position(0))); 
	}

  Vector<double> Sinusoidal::derivative(Vector<double> const & position) const
  { 
    Vector<double> deriv(1);
    deriv(0) = amplitude_*pulsation_*std::sin(pulsation_*position(0));
    return deriv;
  }
  
  double Sinusoidal::laplacian(Vector<double> const & position) const
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
  
  double SumSinusoidal::operator()(Vector<double> const & position) const  { return amplitude_* (std::sin(pulsation_*position(0)) 
    
	+ cos(2 * pulsation_*position(0)) 
	+ cos(3 * pulsation_*position(0)) / 3); }

  Vector<double> SumSinusoidal::derivative(Vector<double> const & position) const
  { 
    Vector<double> deriv(1);
    deriv(0) = amplitude_*pulsation_*(cos(pulsation_*position(0))
	  - 2 * sin(2 * pulsation_*position(0))
	  - sin(3 * pulsation_*position(0)));
    return deriv;
  }

    
  double SumSinusoidal::laplacian(Vector<double> const & position) const
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
  
  double FracSinusoidal::operator()(Vector<double> const & position) const
  { 
    double q = position(0);
    return amplitude_ * cos(2 * M_PI * q) / (2 + sin(2 * M_PI * q));
  }

  Vector<double> FracSinusoidal::derivative(Vector<double> const & position) const
  { 
    double q = position(0);
    Vector<double> deriv(1);
    //deriv(0) = -(2 * M_PI * sin(2 * M_PI * q))/(2 + sin(2 * M_PI * q))
	//	-(2 * M_PI * pow(cos(2 M_PI * q),2))/(sin(2 * M_PI * q)+2)^2;
    deriv(0) = -amplitude_ * (4 * M_PI * sin(2 * M_PI * q) + 2 * M_PI)/pow(sin(2 * M_PI * q) + 2, 2);
    return deriv;
  }

    
  double FracSinusoidal::laplacian(Vector<double> const & position) const
  {
    double q = position(0);    
    return -amplitude_ * 32 * pow(M_PI,2) * pow(sin(M_PI/4-M_PI * q), 3) * sin(M_PI * q + M_PI/4)
	  / pow(sin(2 * M_PI * q) + 2, 3);
  }
  
//#### DoubleWell #####  
  
    DoubleWell::DoubleWell(Input const & input, int const& indexOfReplica):Potential(input, indexOfReplica), height_(input.height()), interWell_(input.interWell())
  {}
  
  double DoubleWell::operator()(Vector<double> const & position) const
  { return height_*pow(position(0)-interWell_/2, 2)*pow(position(0)+interWell_/2, 2); }

  Vector<double> DoubleWell::derivative(Vector<double> const & position) const
  { 
    Vector<double> deriv(1);
    deriv(0) = 4*height_*position(0)*(position(0)-interWell_/2)*(position(0)+interWell_/2);
    return deriv;
  }
  
  //#### HarmonicWell #####  
  
  HarmonicWell::HarmonicWell(Input const & input, int const& indexOfReplica):
    Potential(input, indexOfReplica), 
    stiffness_(input.potentialStiffness())
  {}
  
  
  double HarmonicWell::operator()(Vector<double> const & position) const
  { return stiffness_ / 2 * pow(position(0), 2); }

  Vector<double> HarmonicWell::derivative(Vector<double> const & position) const
  { 
    Vector<double> deriv(1);
    deriv(0) = stiffness_ * position(0);
    return deriv;
  }

  
//#### Harmonic #####  
  
  Harmonic::Harmonic(Input const & input, int const& indexOfReplica):
    Potential(input, indexOfReplica), 
    stiffness_(input.potentialStiffness())
  {}
  
  
  double Harmonic::operator()(Vector<double> const & vecDistance) const
  { return stiffness_ / 2* pow(vecDistance(0) - 1, 2); }

  Vector<double> Harmonic::derivative(Vector<double> const & vecDistance) const
  { 
    Vector<double> deriv(1);
    deriv(0) = stiffness_ * (vecDistance(0) - 1);
    return deriv;
  }
  
  //#### Rotor #####  
  
  Rotor::Rotor(Input const & input, int const& indexOfReplica):
    Potential(input, indexOfReplica)
  {}
  
  
  double Rotor::operator()(Vector<double> const & vecDistance) const
  { 
    //return pow(vecDistance(0) - 1, 2);
    return 1 - cos(vecDistance(0));
  }

  Vector<double> Rotor::derivative(Vector<double> const & vecDistance) const
  { 
    Vector<double> deriv(1);
    deriv(0) = sin(vecDistance(0));
    //deriv(0) = 2 * (1 - vecDistance(0));
    return deriv;
  }

  double Rotor::laplacian(Vector<double> const & vecDistance) const
  { 
    return cos(vecDistance(0));
  }
  
    
//#### Quadratic #####  
  
  Quadratic::Quadratic(Input const & input, int const& indexOfReplica):
    Potential(input, indexOfReplica), 
    stiffness_(input.potentialStiffness()),
    alpha_(input.potentialAlpha()),
    beta_(input.potentialBeta())
  {}
  
  
  double Quadratic::operator()(Vector<double> const & vecDistance) const
  { return stiffness_/2 * pow(vecDistance(0), 2) + alpha_/3 * pow(vecDistance(0), 3) + beta_/4 * pow(vecDistance(0), 4); }

  Vector<double> Quadratic::derivative(Vector<double> const & vecDistance) const
  { 
    Vector<double> deriv(1);
    deriv(0) = stiffness_ * vecDistance(0) + alpha_ * pow(vecDistance(0), 2) + beta_ * pow(vecDistance(0), 3);
    return deriv;
  }
}


#endif
