#include "simol/statphys/potential/AllPotentials.hpp"

namespace simol
{
  
  Potential* createPotential(Input const& input, string potName)
  {
    if (potName == "None")
      return new NoPotential(input);
    else if (potName == "Sinusoidal")
      return new Sinusoidal(input);
    else if (potName == "SumSinusoidal")
      return new SumSinusoidal(input);
    else if (potName == "FracSinusoidal")
      return new FracSinusoidal(input);
    else if (potName == "DoubleWell")
      return new DoubleWell(input);
    else if (potName == "Harmonic")
      return new Harmonic(input);
    else if (potName == "Rotor")
      return new Rotor(input);
    else if (potName == "FPU")
      return new FPU(input);
    else if (potName == "SpaceSinus")
      return new SpaceSinus(input);
    else if (potName == "LennardJones")
      return new LennardJones(input);
    else if (potName == "WCA")
      return new WCA(input);
    else if (potName == "SoftDPD")
      return new SoftDPD(input);
    else if (potName == "TwoTypes")
      return new TwoTypes(input);
    else if (potName == "HarmonicFE")
      return new HarmonicFE(input);
    else if (potName == "DoubleWellFE")
      return new DoubleWellFE(input);
    else
      throw runtime_error(potName +" is not a valid potential !");
  }

  /*Potential* createPotential(Input const& input)
  {
    return createPotential(input, input.potentialName());
  }
  
  ///
  ///Returns a pointer toward a 1D potential, suitable for Galerkin computation
  Potential* createGalerkinPotential(Input const& input)
  {
    if (input.potentialName() == "TwoTypes")
      return createPotential(input, input.secondPotentialName());
    else
      return createPotential(input, input.potentialName());
  }*/
  
  
  TwoTypes::TwoTypes(Input const & input): 
    Potential(input),
    interactionRatio_(input.interactionRatio())
  {
    firstPot_ = createPotential(input, input.firstPotentialName());
    secondPot_= createPotential(input, input.secondPotentialName());
  }
  
  double TwoTypes::operator()(double position, int type) const
  {
    if (type == 0)
      return firstPot_->value(position);
    else if (type == 1)
      return interactionRatio_ * firstPot_->value(position);
    else if (type == 2)
      return secondPot_->value(position);
    else
      throw runtime_error("Potential type int invalid !");
  }

  DVec TwoTypes::gradient(double position, int type) const
  {
    if (type == 0)
      return firstPot_->gradient(position);
    else if (type == 1)
      return interactionRatio_ * firstPot_->gradient(position);
    else if (type == 2)      
      return secondPot_->gradient(position);
    else
      throw runtime_error("Potential type int invalid !");
  }
  
  DVec TwoTypes::gradient(double /*position*/) const
  {
    throw runtime_error("An interaction type should be specified for a call to a TwoTypes potential !");
  }
  
  DVec TwoTypes::potentialForce(double position, int type) const
  {
    return -gradient(position, type);
  }
  
  /*DVec TwoTypes::potentialForce(double position, int type) const
  {
    if (type == 0)
      return firstPot_->potentialForce(position);
    else (type == 1)
      return interactionRatio_ * secondPot_->gradient(position);
    else (type == 2)
      return secondPot_->potentialForce(position);
    else
      throw runtime_error("Potential type int invalid !");
  }
  
  DVec TwoTypes::potentialForce(double position) const
  {
    throw runtime_error("An interaction type should be specified for a call to a TwoTypes potential !");
  }*/
  
  ///
  /// Ratio used in the rejection method
  double TwoTypes::shiftToHarmonic(int type) const
  {
    if (type == 0)
      return firstPot_->shiftToHarmonic();
    else
      return secondPot_->shiftToHarmonic();
  }
  
  double TwoTypes::drawLaw(double localBeta, std::shared_ptr<RNG>& rng, int type) const
  {
    if (type < 2)
      return firstPot_->drawLaw(localBeta, rng, type);
    else
      return secondPot_->drawLaw(localBeta, rng, type);
  }
  
  
  NoPotential::NoPotential(Input const& input):
    Potential(input)
  {}

}
