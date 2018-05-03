#include "simol/statphys/system/Bicolor.hpp"

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
    else if (systemSubtype() == "ColorDrift" || systemSubtype() == "NortonColorDrift")
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        getParticle(iOfParticle).type() = 2*((2*iOfParticle) < nbOfParticles()) - 1;    
    else
      throw runtime_error(systemSubtype() + " is not a valid system subtype !");
  }
  
  void Bicolor::computeAllForces()
  {
    //cout << endl << "------------------computeAllForces----------------------" << endl;
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "1 : i " << iOfParticle << " " << getParticle(iOfParticle).type() << endl;

    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).resetForce(externalPotential());
    
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "2 : i " << iOfParticle << ", p " << getParticle(iOfParticle).position().adjoint() << ", e " << getParticle(iOfParticle).potentialEnergy() << ", f " << getParticle(iOfParticle).force().adjoint() << endl;
    
    
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
    
    /*double alpha3=pairPotential().interactionRatio();
    lagrangeMultiplier() = - getParticle(0).force(0) / alpha3;
    getParticle(0).force(0) *= (2-alpha3)/alpha3;      */
    
  }
  
  void Bicolor::enforceConstraint(double drift, DynamicsParameters const& /*dynaPara*/, bool updateLagrangeMultiplier)
  {
    
    /*double relativeDrift = 0;
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      relativeDrift += getParticle(iOfParticle).type() * getParticle(iOfParticle).velocity(0);
    relativeDrift /= sqrt((double) nbOfParticles());*/
    double relativeDrift = velocity();
    
    //cout << "relativeDrift before = " << relativeDrift << " compared to drift = " << drift << endl;
    
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
      //return sumVelocity;
      sumVelocity /= (double) nbOfParticles();
      //return (sumVelocity + externalPotential_->nonEqAmplitude() /(nbOfParticles()-1.)) * (nbOfParticles()-1.) / nbOfParticles();
      //return (sumVelocity + externalPotential_->nonEqAmplitude() /(dynaPara.gamma() * nbOfParticles()-1.)) * (nbOfParticles()-1.) / nbOfParticles();
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