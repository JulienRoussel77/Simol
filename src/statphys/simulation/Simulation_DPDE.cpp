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
    //-- Verlet part --
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.secondPartThermalization(syst(iOfParticle));
  }

  //------------- DPDE --------------------
  
  void simulate(DPDE& dyna, System& syst)
  {
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletSecondPart(syst(iOfParticle));
    
    //---------- Isolated systems -------------
    if (syst.name() == "Isolated")
      {
	//-- fluctuation/dissipation --
	for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
	  //dyna.energyReinjection(syst(iOfParticle));  // integration of p at fixed gamma + energy reinjection
	  dyna.metropolizedEnergyReinjection(syst(iOfParticle));  // Metropolis correction of effective dynamics on p 
      }
    //------------- NBody --------------------
    else if (syst.name() == "NBody")
      {
	//-- fluctuation/dissipation --
	dyna.fluctuationDissipation(syst);
	//-- thermal conduction --
	//dyna.thermalConduction(syst);
      }
    else 
      throw std::runtime_error("The system is not implemented for DPDE dynamics!");
    
  }

}
  
 
