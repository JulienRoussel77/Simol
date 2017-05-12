#ifndef SIMOL_ALLPOTENTIALS_HPP
#define SIMOL_ALLPOTENTIALS_HPP

#include "DoubleWell.hpp"
#include "FPU.hpp"
#include "FracSinusoidal.hpp"
#include "Harmonic.hpp"
#include "LennardJones.hpp"
#include "Rotor.hpp"
#include "Sinusoidal.hpp"
#include "SpaceSinus.hpp"
#include "SumSinusoidal.hpp"
#include "SoftDPD.hpp"

namespace simol
{
  Potential* createPotential(Input const& input, string potName);
  Potential* createPotential(Input const& input);
  Potential* createGalerkinPotential(Input const& input);
  
  class TwoTypes : public Potential
  {
    public:
      TwoTypes(Input const& input);
      double operator()(double position, int type) const;
      DVec gradient(double position, int type) const;
      DVec gradient(double position) const;
      DVec potentialForce(double position, int type) const;
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
    double operator()(double /*position*/) const {return 0;}
    DVec gradient(double /*position*/) const {return DVec::Zero(dimension());}
    //DVec potentialForce(double /*position*/, int /*type*/) const {return DVec::Zero(1);};
    //DVec potentialForce(double position) const;
    
    virtual double shiftToHarmonic(int /*type*/) const {return 0;};
    virtual double drawLaw(double /*localBeta*/, std::shared_ptr<RNG>& /*rng*/, int /*type*/) const {return 0;}
  };

}

#endif
