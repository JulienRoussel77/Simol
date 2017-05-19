#ifndef POTENTIAL_IMPL_HPP
#define POTENTIAL_IMPL_HPP

#include "simol/statphys/potential/Potential.hpp"

#include <cmath>

using std::cout;
using std::endl;

namespace simol
{

  Potential::Potential(){}

  Potential::Potential(Input const & input):
    nonEqForce_(DVec::Zero(input.dimension())),
    center_(input.potentialCenter()),
    domainSize_(std::numeric_limits<double>::infinity()),
    dimension_(input.dimension())
  {
    nonEqForce_(0) = input.nonEqForce();
    //if (nonEqForce_(0) != 0)
    // cout << "nonEqForce = " << nonEqForce_ << endl;
    cout << "Potential::Potential" << endl;
    cout << "nonEqForce : " << nonEqForce().adjoint() << endl;
    
  }

  ///
  ///Read-only accessor for the external force
  DVec const& Potential::nonEqForce() const
  {
    return nonEqForce_;
  }
  ///
  /// Write-read accessor for the external force
  DVec& Potential::nonEqForce()
  {
    return nonEqForce_;
  }
  ///
  ///Read-only accessor for the i-th component of the external force
  double const& Potential::nonEqForce(int const& i) const
  {
    return nonEqForce_(i);
  }
  ///
  ///Write-read accessor for the i-th component of the external force
  double& Potential::nonEqForce(int const& i)
  {
    return nonEqForce_(i);
  }
  
  double const& Potential::parameter1() const
  {
    throw std::runtime_error("parameter1 not defined for this potential");
  }
  
  double const& Potential::parameter2() const
  {
    throw std::runtime_error("parameter2 not defined for this potential");
  }
  
  const double& Potential::domainSize() const
  {return domainSize_;}
  
  double& Potential::domainSize()
  {return domainSize_;}
      
  const int& Potential::dimension() const
  {return dimension_;}   
      

  double Potential::operator()(DVec const& /*position*/) const
  {
    throw runtime_error("Potential::operator()(DVec const& position) not implemented !");
    //cout << "Potential::operator()(DVec const& position)" << endl;
    //return operator()(position(0));
  }

  double Potential::operator()(double /*position*/) const
  {
    throw runtime_error("Potential::operator()(double position) not implemented !");
    //cout << "Potential::operator()(double position)" << endl;
    //return operator()(DVec::Constant(dimension(), position));
    //return operator()(DVec::Constant(1,1, position));
  }

  double Potential::value(DVec const& position) const
  {
    //cout << "Potential::value(DVec const& position)" << endl;
    return operator()(position);
  }

  double Potential::value(double position) const
  {
    //cout << "Potential::value(double position)" << endl;
    return operator()(position);
  }
  
  double Potential::symmetricValue(double position) const
  {
    return (value(position) + value(2*center_ - position))/2;
  }
  
  double Potential::skewsymmetricValue(double position) const
  {
    return (value(position) - value(2*center_ - position))/2;
  }

  DVec Potential::gradient(DVec const& position) const
  {
    //throw runtime_error("Potential::gradient(DVec const& position) not implemented !");
    //cout << "Potential::gradient(DVec const& position)" << endl;
    return DVec::Constant(1, 1, scalarGradient(position(0)));
  }

  double Potential::scalarGradient(double /*position*/) const
  {
    throw runtime_error("Potential::scalarGradient(double position) not implemented !");
    //cout << "Potential::gradient(double position)" << endl;
    //return gradient(DVec::Constant(dimension(), position));
  }

  DVec Potential::totalForce(DVec const& position) const
  {
    //cout << "Potential::totalForce" << endl;
    //cout << nonEqForce_ << " - " << gradient(position) << " = " << nonEqForce_ - gradient(position) << endl;
    return nonEqForce_ - gradient(position);
  }

  /*double Potential::scalarTotalForce(double position) const
  {
    return nonEqForce_ - gradient(position);
  }*/

  DVec Potential::potentialForce(DVec const& position) const
  {
    return - gradient(position);
  }

  /*DVec Potential::potentialForce(double position) const
  {
    return - gradient(position);
  }*/
  
  double Potential::scalarPotentialForce(double position) const
  {
    return - scalarGradient(position);
  }

  double Potential::laplacian(DVec const& position) const
  {
    return laplacian(position(0));
  }

  double Potential::laplacian(double position) const
  {
    return laplacian(DVec::Constant(dimension(), position));
  }
  
  ///
  /// Ratio used in the rejection method
  double Potential::shiftToHarmonic() const
  {
    throw std::invalid_argument("Potential::shiftToHarmonic : Function undefined");
  }
  
  double Potential::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
  {
    return drawLaw(localBeta, rng, 0);
  }
  
  
  double Potential::operator()(DVec const & position, int /*type*/) const
  {return operator()(position);}  
  double Potential::operator()(double position, int /*type*/) const
  {return operator()(position);}
  double Potential::value(DVec const& position, int /*type*/) const
  {return value(position);}
  double Potential::value(double position, int /*type*/) const
  {return value(position);}
  DVec Potential::gradient(DVec const & position, int /*type*/) const
  {return gradient(position);}
  double Potential::scalarGradient(double position, int /*type*/) const
  {return scalarGradient(position);}
  DVec Potential::totalForce(DVec const & position, int /*type*/) const
  {return totalForce(position);}
  /*DVec Potential::totalForce(double position, int type) const
  {return totalForce(position);}*/
  DVec Potential::potentialForce(DVec const & position, int /*type*/) const
  {return potentialForce(position);}
  double Potential::scalarPotentialForce(double position, int /*type*/) const
  {return scalarPotentialForce(position);}
  double Potential::laplacian(DVec const & position, int /*type*/) const
  {return laplacian(position);}
  double Potential::laplacian(double position, int /*type*/) const
  {return laplacian(position);}
      


  ///
  /// Ratio used in the rejection method
  double Potential::shiftToHarmonic(int /*type*/) const
  {
    return shiftToHarmonic();
  }
  
  /*///
  ///sampling using the rejection method
  double Potential::drawLaw(double localBeta, std::shared_ptr<RNG>& rng, int type) const
  {
    double ratio = shiftToHarmonic(type);
    bool reject = true;
    double xdraw, udraw;
    cout << "Rejection method with a gaussian centered on " << center_ << " with std dev of 1 and ratio " << ratio << endl;
    while (reject)
    {
      xdraw = center_ + rng->scalarGaussian() / sqrt(localBeta);
      cout << "xdraw : " << xdraw << endl;
      
      udraw = rng->scalarUniform();

      reject = (udraw > exp(- localBeta * (value(xdraw, type) - pow(xdraw-center_, 2) / 2 + ratio)));
      cout << udraw << " compared to " << exp(- localBeta * (value(xdraw, type) - pow(xdraw-center_, 2) / 2 + ratio)) << endl;
      assert(exp(-localBeta * (pow(xdraw-center_, 2) / 2 - ratio)) >= exp(- localBeta * value(xdraw)));
    }
    return xdraw;
  }*/
  
  ///
  ///sampling using the inverse tranform method
  double Potential::drawLaw(double localBeta, std::shared_ptr<RNG>& rng, int type) const
  {
    double step = 1e-4;
    double qMin = -10;
    int nbOfNodes = (int)(-2*qMin /step);
    double partitionFunction = 0;
    for (int iOfNode = 0; iOfNode < nbOfNodes; iOfNode++)
    {
      double q = qMin + iOfNode * step;
      partitionFunction += exp(-localBeta * value(q, type));
      //cout << iOfNode << " " << q << " " << exp(-localBeta * value(q, type)) << " " << partitionFunction << endl;
    }
    partitionFunction *= step;
    double randVal = rng->scalarUniform();
    cout << "Inverse transform method : U = " << randVal << endl;
    double objective = randVal * partitionFunction / step;
    double cumulant = 0;
    for (int iOfNode = 0; iOfNode < nbOfNodes; iOfNode++)
    {
      double q = qMin + iOfNode * step;
      cumulant += exp(-localBeta * value(q, type));
      if (cumulant >= objective)
        return q + step/2;
    }
    return qMin + nbOfNodes * step;
  }
  
  double Potential::harmonicForce(double /*dist*/) const
  {
    throw std::invalid_argument("Potential::harmonicForce : Function undefined");
  }
  
  double Potential::harmonicStiffness() const
  {
    throw std::invalid_argument("Potential::harmonicStiffness : Function undefined");
  }

  double Potential::harmonicEquilibrium() const
  {
    throw std::invalid_argument("Potential::harmonicEquilibrium: Function undefined");
  }
  
  double Potential::harmonicFrequency() const
  {
    throw std::invalid_argument("Potential::harmonicFrequency: Function undefined");
  }
  
  DVec Potential::polyCoeffs() const
  {
    throw std::invalid_argument("Potential::polynomialCoeffs: Function undefined");
  }
  
  DVec Potential::polyDCoeffs() const
  {
    return polynomialDerivative(polyCoeffs());
  }
  
  DVec Potential::polyDDCoeffs() const
  {
    return polynomialDerivative(polyDDCoeffs());
  }
  
  ///
  /// In the case when V(q) = a/q + poly(q) returns "a"
  double Potential::inverseCoeff() const
  {
    return 0;
  }
}
#endif
