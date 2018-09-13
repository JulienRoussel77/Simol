#ifndef SIMOL_POTENTIAL_HPP
#define SIMOL_POTENTIAL_HPP

#include "simol/Tools.hpp"
//#include "simol/core/linalg/Vector.hpp"
#include "simol/input/Input.hpp"
#include "simol/core/random/RNG.hpp"


namespace simol
{

  class Potential
  {

    public:
      virtual ~Potential() {};
      Potential(Input const& input);
      virtual string classname() const {return "Potential";}
      double const& nonEqAmplitude() const;
      DVec& nonEqForce() ;
      DVec const& nonEqForce() const;
      double& nonEqForce(const int& i);
      double const& nonEqForce(const int& i) const;
      virtual double const& parameter1() const;
      virtual double const& parameter2() const;
      virtual const double& domainSize() const;
      virtual double& domainSize();
      virtual int const& dimension() const;
      virtual bool const& flag() const;
      virtual double const& interactionRatio() const;
      virtual double const& cutOffRadius() const;
      
      virtual double operator()(DVec const & position) const;
      virtual double operator()(double position) const;
      virtual double value(DVec const& position) const;
      virtual double value(double position) const;
      virtual DVec gradient(DVec const & position) const;
      virtual double scalarGradient(double position) const;
      virtual DVec totalForce(DVec const & position) const;
      virtual DVec potentialForce(DVec const & position) const;
      virtual double scalarPotentialForce(double position) const;
      virtual double laplacian(DVec const & position) const;
      virtual double laplacian(double position) const;
      
      virtual double shiftToHarmonic() const;
      virtual double drawLaw(double /*localBeta*/, std::shared_ptr<RNG>& /*rng*/) const;
 
      virtual double operator()(DVec const & position, int type) const;
      virtual double operator()(double position, int type) const;
      virtual double value(DVec const& position, int type) const;
      virtual double value(double position, int type) const;
      virtual double symmetricValue(double position) const;
      virtual double skewsymmetricValue(double position) const;
      virtual DVec gradient(DVec const & position, int type) const;
      virtual double scalarGradient(double position, int type) const;
      virtual DVec totalForce(DVec const & position, int type) const;
      //virtual DVec totalForce(double position, int type) const;
      virtual DVec potentialForce(DVec const & position, int type) const;
      virtual double scalarPotentialForce(double position, int type) const;
      virtual double laplacian(DVec const & position, int type) const;
      virtual double laplacian(double position, int type) const;
      
      virtual double shiftToHarmonic(int type) const;
      virtual double drawLaw(double /*localBeta*/, std::shared_ptr<RNG>& /*rng*/, int type) const;
      
      virtual double harmonicForce(double dist) const;
      virtual double harmonicStiffness() const;
      virtual double harmonicEquilibrium() const;
      virtual double harmonicFrequency() const;
      
      virtual DVec polyCoeffs() const;
      virtual DVec polyDCoeffs() const;
      virtual DVec polyDDCoeffs() const;
      virtual double inverseCoeff() const;
      
      virtual double marginalWithoutCoupling(DVec const& position) const;
      virtual double gradientCoupling(DVec const& position) const;
    protected:
      Potential();
      // amplitude eta of the driving force
      double nonEqAmplitude_;
      // direction of the driving force
      DVec nonEqForce_;
      // center of the potential, allows to shift the system (useful for the colloid interaction)
      double center_;
      // useful for example for the cosine 1D potential
      double domainSize_;
      // dimension of the phase space
      int dimension_;
      // what is this ?
      bool flag_;
      // allows to switch off a part of the potential (eg for the colloid)
      double interactionRatio_;
      // maximum range of the pair interaction
      double cutOffRadius_;
  };

}
#endif
