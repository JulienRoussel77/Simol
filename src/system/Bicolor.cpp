#include "simol/system/Bicolor.hpp"

namespace simol
{
  //### Bicolor ###
  
  ///
  /// type = -1 : drift to the left
  /// type = 0  : no drift
  /// type =  1 : drift to the right
  Bicolor::Bicolor(Input const& input):
    NBody(input)
  {
    if (systemSubtype() == "None")
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        getParticle(iOfParticle).type() = 0; 
    else if (systemSubtype() == "OneDrift" || systemSubtype() == "NortonOneDrift")
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        getParticle(iOfParticle).type() = (iOfParticle == 0); 
    else if (systemSubtype() == "TwoDrift")
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        getParticle(iOfParticle).type() = (iOfParticle < 2); 
    // Beware, ColorDrift only works for an even number of particles
    else if (systemSubtype() == "ColorDrift" || systemSubtype() == "NortonColorDrift")
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        getParticle(iOfParticle).type() = 2*((2*iOfParticle) < nbOfParticles()) - 1;    
    else
      throw runtime_error(systemSubtype() + " is not a valid system subtype !");
  }
  
  void Bicolor::computeAllForces()
  {
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).resetForce(externalPotential());
    
    if (doCells_)
    {
      //-- reinitialize cells before looping on the pair interactions --
      reinitializeCells();
      //-- compute the interactions --
      for (ParticlePairIterator it = pairBegin(); !pairFinished(it); incrementePairIterator(it))
        interaction(it.particle1(), it.particle2());
    }
    else
    {
      //-- no cell method: std double loop --
      for (int i = 0; i < nbOfParticles(); i++)
        for (int j = i + 1; j < nbOfParticles(); j++)
          interaction(getParticle(i), getParticle(j));
    }        
  }
  
  ///
  /// Projects the momenta according to the Norton algorithm so that velocity() == drift
  void Bicolor::enforceConstraint(double drift, DynamicsParameters const& /*dynaPara*/, bool updateLagrangeMultiplier)
  {
    // we call velocity() which returns the linear combination of velocities which is constrained
    double relativeDrift = velocity();
    
    double localLagrangeMultiplier = getParticle(0).mass() * (drift - relativeDrift);

    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).momentum(0) += localLagrangeMultiplier * getParticle(iOfParticle).type();
    
    if (updateLagrangeMultiplier)
      lagrangeMultiplier() += localLagrangeMultiplier;
  }

  
  ///
  ///Computes the instant value of the observable (averaged) velocity F^\top p/(m*n_F)
  double Bicolor::velocity() const
  {
    if (systemSubtype() == "OneDrift" || systemSubtype() == "None")
      return getParticle(0).velocity(0);
    else if (systemSubtype() == "TwoDrift")
      return (getParticle(0).velocity(0) + getParticle(1).velocity(0))/2.;
    else if (systemSubtype() == "ColorDrift" || systemSubtype() == "NortonColorDrift")
    {
      double sumVelocity = 0;
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        sumVelocity += getParticle(iOfParticle).type() * getParticle(iOfParticle).velocity(0);
      
      sumVelocity /= (double) nbOfParticles();
      
      return sumVelocity;
    }
    else
      throw runtime_error(systemSubtype() + " is not a valid system subtype in Bicolor::velocity() !");
  }
  
  ///
  ///Computes the instant value of the observable (averaged) force F^\top \nabla V
  double Bicolor::force() const
  {
    if (systemSubtype() == "OneDrift" || systemSubtype() == "None")
      return getParticle(0).force(0);
    else if (systemSubtype() == "TwoDrift")
      return getParticle(0).force(0) + getParticle(1).force(0);
    else if (systemSubtype() == "ColorDrift" || systemSubtype() == "NortonColorDrift")
    {
      double sumForces = 0;
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        sumForces += getParticle(iOfParticle).type() * getParticle(iOfParticle).force(0);
      return sumForces;
    }
    else
      throw runtime_error(systemSubtype() + " is not a valid system subtype in Bicolor::force() !");
  }
  
}