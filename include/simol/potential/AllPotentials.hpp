#ifndef SIMOL_ALLPOTENTIALS_HPP
#define SIMOL_ALLPOTENTIALS_HPP

#include "DoubleWell.hpp"
#include "FPU.hpp"
#include "FracSinusoidal.hpp"
#include "Harmonic.hpp"
#include "LennardJones.hpp"
#include "Rotor.hpp"
#include "Sinusoidal.hpp"
#include "SpaceSine.hpp"
#include "Shearing.hpp"
#include "SumSinusoidal.hpp"
#include "SoftDPD.hpp"
#include "Drift.hpp"
#include "Coulomb.hpp"

namespace simol
{
  Potential* createPotential(Input const& input, string potName);
  //Potential* createPotential(Input const& input);
  //Potential* createGalerkinPotential(Input const& input);
  
  class TwoTypes : public Potential
  {
    public:
      TwoTypes(Input const& input);
      virtual string classname() const {return "TwoTypes";}
      double operator()(double position, int type) const;
      double scalarGradient(double position, int type) const;
      double scalarGradient(double position) const;
      double scalarPotentialForce(double position, int type) const;
      //DVec potentialForce(double position) const;
      
      virtual double shiftToHarmonic(int type) const;
      virtual double drawLaw(double localBeta, std::shared_ptr<RNG>& rng, int type) const;
    private:
      Potential* firstPot_;
      Potential* secondPot_;
      double interactionRatio_;
  };
  
  class NoPotential : public Potential
  {
  public:
    NoPotential(Input const& input);
    virtual string classname() const {return "NoPotential";}
    virtual double operator()(DVec const& /*position*/) const;
    virtual double operator()(double /*position*/) const;
    virtual DVec gradient(DVec const& /*position*/) const;
    virtual double scalarGradient(double /*position*/) const;

    //DVec potentialForce(double /*position*/, int /*type*/) const {return DVec::Zero(1);};
    //DVec potentialForce(double position) const;
    
    virtual double shiftToHarmonic(int /*type*/) const {return 0;};
    virtual double drawLaw(double /*localBeta*/, std::shared_ptr<RNG>& /*rng*/, int /*type*/) const {return 0;}
  };

}

#endif
