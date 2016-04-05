#ifndef POTENTIAL_IMPL_HPP
#define POTENTIAL_IMPL_HPP

#include "Potential.hpp"

#include <cmath>

using std::cout; 
using std::endl; 

namespace simol
{

  Potential::Potential(){}
  
  Potential::Potential(Input const & /*input*/){}
  
  Potential* createPotential(Input const& input)
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
    else if (input.potentialName() == "SpaceSinus")
      return new SpaceSinus(input);
    else if (input.potentialName() == "LennardJones")
      return new LennardJones(input);
    else
      std::cout << input.potentialName() << " is not a valid potential !" << std::endl;
    
    return 0;
  }
  
  double Potential::operator()(Vector<double> const& position) const
  {
    //cout << "Potential::operator()(Vector<double> const& position)" << endl;
    return operator()(position(0));
  }
  
  double Potential::operator()(double position) const
  {
    //cout << "Potential::operator()(double position)" << endl;
    return operator()(Vector<double>(1, position));
  }
  
  double Potential::value(Vector<double> const& position) const
  {
    //cout << "Potential::value(Vector<double> const& position)" << endl;
    return operator()(position);
  }
  
  double Potential::value(double position) const
  {
    //cout << "Potential::value(double position)" << endl;
    return operator()(position);
  }
  
  Vector<double> Potential::derivative(Vector<double> const& position) const
  { 
    return derivative(position(0)); 
  }
  
  Vector<double> Potential::derivative(double position) const
  { 
    return -derivative(Vector<double>(1, position)); 
  }
  
  Vector<double> Potential::force(Vector<double> const& position) const
  { 
    return -derivative(position); 
  }
  
  Vector<double> Potential::force(double position) const
  { 
    return -derivative(position); 
  }
  
  double Potential::laplacian(Vector<double> const& position) const
  { 
    return laplacian(position(0)); 
  }
  
  double Potential::laplacian(double position) const
  { 
    return laplacian(Vector<double>(1, position));
  }
  
  double Potential::ratioToHarmonic() const 
  {throw std::invalid_argument("Potential::ratioToHarmonic : Function undefined");}
  
  double Potential::drawLaw(double /*localBeta*/, std::shared_ptr<RNG>& /*rng*/) const
  {throw std::invalid_argument("Potential::drawLaw : Function undefined");}
  
  //#### Sinusoidal #####
  
  Sinusoidal::Sinusoidal(Input const & input):
    Potential(input), 
    amplitude_(input.amplitude()), 
    pulsation_(2*M_PI/input.length())
  {}
  
  double Sinusoidal::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
  {
    bool reject = true;
    double xdraw, udraw;
    int count = 0;
    while (reject)
      {
	xdraw = -M_PI + 2 * rng->scalarUniform() * M_PI;
	udraw = rng->scalarUniform();
	
	reject = (udraw > exp(- localBeta * (value(xdraw) + 2)));
	count++;
      }
    return xdraw;
  }
  
  double Sinusoidal::operator()(double position) const
  { 
    return amplitude_* (1-cos(pulsation_*position)); 
  }
  
  Vector<double> Sinusoidal::derivative(double position) const
  { 
    return Vector<double>(1, amplitude_*pulsation_*sin(pulsation_*position));
  }
  
  double Sinusoidal::laplacian(double position) const
  {
    return amplitude_ * pow(pulsation_, 2) * cos(pulsation_ * position);
  }
  
    
//#### SumSinusoidal #####
  
  SumSinusoidal::SumSinusoidal(Input const & input):
    Potential(input), 
    amplitude_(input.amplitude()), 
    pulsation_(2*M_PI/input.length())
  {}
  
  double SumSinusoidal::operator()(double position) const  
  { 
		return amplitude_* (std::sin(pulsation_*position) 
			+ cos(2 * pulsation_*position) 
			+ cos(3 * pulsation_*position) / 3); 
	}

  Vector<double> SumSinusoidal::derivative(double position) const
  { 
    Vector<double> deriv(1);
    deriv(0) = amplitude_*pulsation_*(cos(pulsation_*position)
	  - 2 * sin(2 * pulsation_*position)
	  - sin(3 * pulsation_*position));
    return deriv;
  }

    
  double SumSinusoidal::laplacian(double position) const
  {
    return amplitude_ * pow(pulsation_, 2) * (-sin(pulsation_ * position)
	      - 4 * cos(2 * pulsation_ * position)
	      - 3 * cos(3 * pulsation_ * position));
  }
  
  
  //#### FracSinusoidal #####
  
  FracSinusoidal::FracSinusoidal(Input const & input):
    Potential(input),
    amplitude_(input.amplitude()),
    pulsation_(2*M_PI/input.length())
  {}
  
  double FracSinusoidal::operator()(double position) const
  { 
    return amplitude_ * cos(2 * M_PI * position) / (2 + sin(2 * M_PI * position));
  }

  Vector<double> FracSinusoidal::derivative(double position) const
  { 
    Vector<double> deriv(1);
    //deriv(0) = -(2 * M_PI * sin(2 * M_PI * position))/(2 + sin(2 * M_PI * position))
	//	-(2 * M_PI * pow(cos(2 M_PI * position),2))/(sin(2 * M_PI * position)+2)^2;
    deriv(0) = -amplitude_ * (4 * M_PI * sin(2 * M_PI * position) + 2 * M_PI)/pow(sin(2 * M_PI * position) + 2, 2);
    return deriv;
  }

    
  double FracSinusoidal::laplacian(double position) const
  {
    return -amplitude_ * 32 * pow(M_PI,2) * pow(sin(M_PI/4-M_PI * position), 3) * sin(M_PI * position + M_PI/4)
	  / pow(sin(2 * M_PI * position) + 2, 3);
  }
  
//#### DoubleWell #####  
  
    DoubleWell::DoubleWell(Input const & input):Potential(input), height_(input.height()), interWell_(input.interWell())
  {}
  
  double DoubleWell::operator()(double position) const
  { return height_*pow(position-interWell_/2, 2)*pow(position+interWell_/2, 2); }

  Vector<double> DoubleWell::derivative(double position) const
  { 
    Vector<double> deriv(1);
    deriv(0) = 4*height_*position*(position-interWell_/2)*(position+interWell_/2);
    return deriv;
  }
  
  //#### HarmonicWell #####  
  
  HarmonicWell::HarmonicWell(Input const & input):
    Potential(input), 
    stiffness_(input.potentialStiffness())
  {}
  
  
  double HarmonicWell::operator()(double position) const
  { return stiffness_ / 2 * pow(position, 2); }

  Vector<double> HarmonicWell::derivative(double position) const
  { 
    Vector<double> deriv(1);
    deriv(0) = stiffness_ * position;
    return deriv;
  }

  
//#### Harmonic #####  
  
  Harmonic::Harmonic(Input const & input):
    Potential(input), 
    stiffness_(input.potentialStiffness())
  {}
  
  
  double Harmonic::operator()(double distance) const
  { return stiffness_ / 2* pow(distance - 1, 2); }

  Vector<double> Harmonic::derivative(double distance) const
  { 
    Vector<double> deriv(1);
    deriv(0) = stiffness_ * (distance - 1);
    return deriv;
  }
  
  double Harmonic::laplacian(double /*distance*/) const
  { 
    return stiffness_;
  }
  
  double Harmonic::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
	{
		return rng->scalarGaussian() / sqrt(localBeta);
	}
  
  //#### Rotor #####  
  
  Rotor::Rotor(Input const & input):
    Potential(input)
  {}
  
  
  double Rotor::operator()(double distance) const
  { 
    //return pow(distance - 1, 2);
    return 1 - cos(distance);
  }

  Vector<double> Rotor::derivative(double distance) const
  { 
    Vector<double> deriv(1);
    deriv(0) = sin(distance);
    //deriv(0) = 2 * (1 - distance);
    return deriv;
  }

  double Rotor::laplacian(double distance) const
  {
    return cos(distance);
  }
  
  double Rotor::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
	{
		bool reject = true;
		double xdraw, udraw;
		int count = 0;
		while (reject)
		{
			xdraw = -M_PI + 2 * rng->scalarUniform() * M_PI;
			udraw = rng->scalarUniform();
			
			reject = (udraw > exp(- localBeta * (value(xdraw) + 2)));
			count++;
		}
		return xdraw;
	}
  
    
//#### Quadratic #####  
  
  Quadratic::Quadratic(Input const & input):
    Potential(input), 
    stiffness_(input.potentialStiffness()),
    alpha_(input.potentialAlpha()),
    beta_(input.potentialBeta())
  {
		assert(beta_ > 0 || (beta_ == 0 && alpha_ == 0));
	}
  
  
  double Quadratic::operator()(double distance) const
  { return stiffness_/2 * pow(distance, 2) + alpha_/3 * pow(distance, 3) + beta_/4 * pow(distance, 4); }

  Vector<double> Quadratic::derivative(double distance) const
  { 
    Vector<double> deriv(1);
    deriv(0) = stiffness_ * distance + alpha_ * pow(distance, 2) + beta_ * pow(distance, 3);
    return deriv;
  }
  
  double Quadratic::laplacian(double distance) const
  { 
    return stiffness_ + 2 * alpha_ * distance + 3 * beta_ * pow(distance, 2);
  }
  
  double Quadratic::ratioToHarmonic() const
  {
		assert(stiffness_ == 1);
		if (beta_ > 0)
			return - pow(alpha_, 4) / (12 * pow(beta_, 3));
		else
			return 0;
	}
	
	double Quadratic::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
	{
		double ratio = ratioToHarmonic();
		bool reject = true;
		double xdraw, udraw;
		while (reject)
		{
			xdraw = rng->scalarGaussian() / sqrt(localBeta);
			//cout << ratio << " " << exp(-localBeta * (pow(xdraw, 2)/2 + ratio)) << " " << exp(- localBeta * potential_->value(xdraw)) << endl;
			//cout << xdraw << " " << localBeta * pow(xdraw, 2)/2 - ratio << " >= " << localBeta * potential_->value(xdraw) << endl;
			udraw = rng->scalarUniform();
			
			reject = (udraw > exp(- localBeta * (value(xdraw) + pow(xdraw, 2)/2 + ratio)));
			//cout << reject << " " << xdraw << " " << ydraw << endl << endl;
			assert(exp(-localBeta * (pow(xdraw, 2)/2 + ratio)) >= exp(- localBeta * value(xdraw)));
		}
		return xdraw;
	}
	
	//#### SpaceSinus #####
  
  SpaceSinus::SpaceSinus(Input const & input):
		Potential(input), 
		amplitude_(input.amplitude()), 
		pulsation_(2*M_PI/input.length())
  {}
  
  double SpaceSinus::operator()(Vector<double> const& position) const
  { 
		return amplitude_* (1-cos(pulsation_*position(0))*cos(pulsation_*position(1))*cos(pulsation_*position(2))); 
	}

  Vector<double> SpaceSinus::derivative(Vector<double> const& position) const
  { 
    Vector<double> grad(3);
    grad(0) = sin(pulsation_*position(0))*cos(pulsation_*position(1))*cos(pulsation_*position(2));
    grad(1) = cos(pulsation_*position(0))*sin(pulsation_*position(1))*cos(pulsation_*position(2));
    grad(2) = cos(pulsation_*position(0))*cos(pulsation_*position(1))*sin(pulsation_*position(2));
    grad *= amplitude_*pulsation_;
    return grad;
  }
  
  double SpaceSinus::laplacian(Vector<double> const& position) const
  {
    return 3 * pow(pulsation_, 2) * value(position);
  }
	
  
  //------------- Lennard Jones ---------------------
  LennardJones::LennardJones(Input const & input):
    Potential(input), 
    epsLJ_(input.epsLJ()), 
    sigmaLJ_(input.sigmaLJ())
  {}
  
  double LennardJones::operator()(double const & dist) const
  { 
    return 4 * epsLJ_* (pow(sigmaLJ_/dist,12)-pow(sigmaLJ_/dist,6));
  }
  
  double LennardJones::derivative(double const & dist) const
  { 
    return -24 * epsLJ_* sigmaLJ_ * (2*pow(dist/sigmaLJ_,13)-pow(dist/sigmaLJ_,7));
  }
  
}
#endif
