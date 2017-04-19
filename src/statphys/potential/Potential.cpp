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
    center_(input.potentialCenter())
  {
    nonEqForce_(0) = input.nonEqForce();
    //if (nonEqForce_(0) != 0)
    // cout << "nonEqForce = " << nonEqForce_ << endl;
    
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

  double Potential::operator()(DVec const& position) const
  {
    //cout << "Potential::operator()(DVec const& position)" << endl;
    return operator()(position(0));
  }

  double Potential::operator()(double position) const
  {
    //cout << "Potential::operator()(double position)" << endl;
    return operator()(DVec::Constant(1,1,position));
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
    return gradient(position(0));
  }

  DVec Potential::gradient(double position) const
  {
    return -gradient(DVec::Constant(1,1, position));
  }

  DVec Potential::totalForce(DVec const& position) const
  {
    return nonEqForce_ - gradient(position);
  }

  DVec Potential::totalForce(double position) const
  {
    return nonEqForce_ - gradient(position);
  }

  DVec Potential::potentialForce(DVec const& position) const
  {
    return - gradient(position);
  }

  DVec Potential::potentialForce(double position) const
  {
    return - gradient(position);
  }

  double Potential::laplacian(DVec const& position) const
  {
    return laplacian(position(0));
  }

  double Potential::laplacian(double position) const
  {
    return laplacian(DVec::Constant(1,1, position));
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
  DVec Potential::gradient(double position, int /*type*/) const
  {return gradient(position);}
  DVec Potential::totalForce(DVec const & position, int /*type*/) const
  {return totalForce(position);}
  DVec Potential::totalForce(double position, int /*type*/) const
  {return totalForce(position);}
  DVec Potential::potentialForce(DVec const & position, int /*type*/) const
  {return potentialForce(position);}
  DVec Potential::potentialForce(double position, int /*type*/) const
  {return potentialForce(position);}
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
    double objective = rng->scalarUniform() * partitionFunction / step;
    double cumulant = 0;
    for (int iOfNode = 0; iOfNode < nbOfNodes; iOfNode++)
    {
      double q = qMin + iOfNode * step;
      cumulant += exp(-localBeta * value(q, type));
      if (cumulant >= objective)
        return q + step/2;
    }
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
