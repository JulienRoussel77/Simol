#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{

  void sampleInternalEnergies(DPDE const&, System& syst)
  {
    cout << " - Sampling internal energies..." << endl;
    for (int i = 0; i < syst.nbOfParticles(); i++)
      syst.getParticle(i).internalEnergy() = 1;
  }

  //--- DPDE dynamics ---
  void thermalize(DPDE& dyna, System& syst)
  {
    //-- Verlet part: possibly MTS here as well, as in the function "simulate"--
    for (int k = 0; k < dyna.MTSfrequency(); k++)
      {
	for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
	  dyna.verletFirstPart(syst(iOfParticle));
	syst.computeAllForces();
	for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
	  dyna.verletSecondPart(syst(iOfParticle));
      }
    //-- for the FD part, the timestep is renormalized as the effective timestep MTSfrequency()*timeStep_, as in the function "simulate" --
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.Thermalization(syst(iOfParticle));
    }

}
  
 
