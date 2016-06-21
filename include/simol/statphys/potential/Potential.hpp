#ifndef SIMOL_POTENTIAL_HPP
#define SIMOL_POTENTIAL_HPP

#include "simol/statphys/Tools.hpp"
#include "simol/core/linalg/Vector.hpp"
#include "simol/statphys/input/Input.hpp"
#include "simol/core/random/RNG.hpp"


namespace simol
{

  class Potential
  {

    public:
      virtual ~Potential() {};
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
      virtual double shiftToHarmonic() const;
      virtual double drawLaw(double /*localBeta*/, std::shared_ptr<RNG>& /*rng*/) const;

    protected:
      Potential();
      Vector<double> externalForce_;

  };

}
#endif
