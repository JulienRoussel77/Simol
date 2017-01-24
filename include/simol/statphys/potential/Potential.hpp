#ifndef SIMOL_POTENTIAL_HPP
#define SIMOL_POTENTIAL_HPP

#include "simol/statphys/Tools.hpp"
//#include "simol/core/linalg/Vector.hpp"
#include "simol/statphys/input/Input.hpp"
#include "simol/core/random/RNG.hpp"


namespace simol
{

  class Potential
  {

    public:
      virtual ~Potential() {};
      Potential(Input const& input);
      DVec& externalForce() ;
      DVec const& externalForce() const;
      double& externalForce(const int& i);
      double const& externalForce(const int& i) const;
      virtual double const& parameter1() const;
      virtual double const& parameter2() const;
 
      virtual double operator()(DVec const & position) const;
      virtual double operator()(double position) const;
      double value(DVec const& position) const;
      double value(double position) const;
      virtual DVec gradient(DVec const & position) const;
      virtual DVec gradient(double position) const;
      virtual DVec totalForce(DVec const & position) const;
      virtual DVec totalForce(double position) const;
      virtual DVec potentialForce(DVec const & position) const;
      virtual DVec potentialForce(double position) const;
      virtual double laplacian(DVec const & position) const;
      virtual double laplacian(double position) const;
      virtual double shiftToHarmonic() const;
      virtual double drawLaw(double /*localBeta*/, std::shared_ptr<RNG>& /*rng*/) const;
      virtual double harmonicForce(double dist) const;
      virtual double harmonicStiffness() const;
      virtual double harmonicEquilibrium() const;
      virtual double harmonicFrequency() const;
    protected:
      Potential();
      DVec externalForce_;

  };

}
#endif
