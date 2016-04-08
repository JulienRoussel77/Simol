#ifndef SIMOL_POTENTIAL_HPP
#define SIMOL_POTENTIAL_HPP

#include "Tools.hpp"
#include "core/linalg/Vector.hpp"
#include "Input.hpp"
#include "core/random/RNG.hpp"


namespace simol
{
  //class Potential;

  //Potential* createPotential(Input const& input);
  
  class Potential
  {
    //friend Potential* createPotential(Input const& input);
  protected:
    Potential();
    Vector<double> externalForce_;
  public:
    virtual ~Potential(){};
    Potential(Input const& input);
    Vector<double>& externalForce() ;
    Vector<double> const& externalForce() const;
    double& externalForce(const int& i);
    double const& externalForce(const int& i) const;
    
    virtual double operator()(Vector<double> const & position) const;
    virtual double operator()(double position) const;
    double value(Vector<double> const& position) const;
    double value(double position) const;
    virtual Vector<double> gradient(Vector<double> const & position) const;
    virtual Vector<double> gradient(double position) const;
    virtual Vector<double> totalForce(Vector<double> const & position) const;
    virtual Vector<double> totalForce(double position) const;
    virtual Vector<double> potentialForce(Vector<double> const & position) const;
    virtual Vector<double> potentialForce(double position) const;
    virtual double laplacian(Vector<double> const & position) const;
    virtual double laplacian(double position) const;
    virtual double ratioToHarmonic() const;
    virtual double drawLaw(double /*localBeta*/, std::shared_ptr<RNG>& /*rng*/) const;
    
  };
  
  
  


}
#endif
