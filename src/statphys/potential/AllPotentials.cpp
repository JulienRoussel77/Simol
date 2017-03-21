#include "simol/statphys/potential/AllPotentials.hpp"

namespace simol
{
  
  Potential* createPotential(Input const& input, string potName)
  {
    if (potName == "Sinusoidal")
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
    else if (potName == "SoftDPD")
      return new SoftDPD(input);
    else if (potName == "TwoTypes")
      return new TwoTypes(input);
    else
      throw runtime_error(potName +" is not a valid potential !");

    return 0;
  }

  Potential* createPotential(Input const& input)
  {
    return createPotential(input, input.potentialName());
  }
  
  
  TwoTypes::TwoTypes(Input const & input): 
    Potential(input)
  {
    firstPot_ = createPotential(input, input.firstPotentialName());
    secondPot_= createPotential(input, input.secondPotentialName());
  }
  
  double TwoTypes::operator()(double position, int type) const
  {
    if (type == 0)
      return firstPot_->value(position);
    else
      return secondPot_->value(position);
  }

  DVec TwoTypes::gradient(double position, int type) const
  {
    if (type == 0)
      return firstPot_->gradient(position);
    else
      return secondPot_->gradient(position);
  }

}
