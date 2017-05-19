#include "simol/statphys/system/Bicolor.hpp"

namespace simol
{
  //### Bicolor ###
  
  ///
  /// type = -1 : drift to the left
  /// type = 0  : no drift
  /// type =  1 : drift to the right
  Bicolor::Bicolor(Input const& input):
    NBody(input),
    fixedVelocity_(input.fixedVelocity())
  {
    if (systemSubtype() == "None")
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        getParticle(iOfParticle).type() = 0; 
    else if (systemSubtype() == "OneDrift")
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        getParticle(iOfParticle).type() = (iOfParticle == 0); 
    else if (systemSubtype() == "TwoDrifts")
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        getParticle(iOfParticle).type() = (iOfParticle < 2); 
    else if (systemSubtype() == "ColorDrift")
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        getParticle(iOfParticle).type() = 2*(2*iOfParticle < nbOfParticles()) - 1;    
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
    
    
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "3 : i " << iOfParticle << ", p " << getParticle(iOfParticle).position().adjoint() << ", e " << getParticle(iOfParticle).potentialEnergy() << ", f " << getParticle(iOfParticle).force().adjoint() << endl;

    
  }
  
  ///
  ///Computes the instant value of the observable (averaged) velocity F^\top p/m
  double Bicolor::velocity() const
  {
    if (systemSubtype() == "OneDrift" || systemSubtype() == "None")
      return getParticle(0).velocity(0);
    else if (systemSubtype() == "TwoDrifts")
      return getParticle(0).velocity(0) + getParticle(1).velocity(0);
    else if (systemSubtype() == "ColorDrift")
    {
      double sumVelocity = 0;
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        sumVelocity += getParticle(iOfParticle).type() * getParticle(iOfParticle).velocity(0);
      return sumVelocity;
    }
    else
      throw runtime_error(systemSubtype() + " is not a valid system subtype in Bicolor::velocity() !");
  }
  
}