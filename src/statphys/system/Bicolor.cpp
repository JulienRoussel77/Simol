#include "simol/statphys/system/Bicolor.hpp"

namespace simol
{
  //### Bicolor ###
  
  Bicolor::Bicolor(Input const& input):
    NBody(input),
    fixedVelocity_(input.fixedVelocity())
  {
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).type() = (2*iOfParticle < nbOfParticles());    
  }
  
  void Bicolor::computeAllForces()
  {
    //cout << endl << "------------------computeAllForces----------------------" << endl;
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "1 : i " << iOfParticle << " " << getParticle(iOfParticle).type() << endl;

    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).resetForce(pairPotential());
    
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "2 : i " << iOfParticle << ", p " << getParticle(iOfParticle).position().adjoint() << ", e " << getParticle(iOfParticle).potentialEnergy() << ", f " << getParticle(iOfParticle).force().adjoint() << endl;
    
    
    for (ParticlePairIterator it = pairBegin(); !pairFinished(it); incrementePairIterator(it))
      interaction(it.particle1(), it.particle2());
    
    
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "3 : i " << iOfParticle << ", p " << getParticle(iOfParticle).position().adjoint() << ", e " << getParticle(iOfParticle).potentialEnergy() << ", f " << getParticle(iOfParticle).force().adjoint() << endl;

    
  }
  
}