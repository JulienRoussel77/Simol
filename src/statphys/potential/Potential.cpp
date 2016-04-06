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
  
  Vector<double> Potential::gradient(Vector<double> const& position) const
  { 
    return gradient(position(0)); 
  }
  
  Vector<double> Potential::gradient(double position) const
  { 
    return -gradient(Vector<double>(1, position)); 
  }
  
  Vector<double> Potential::force(Vector<double> const& position) const
  { 
    return -gradient(position); 
  }
  
  Vector<double> Potential::force(double position) const
  { 
    return -gradient(position); 
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

  
}
#endif
