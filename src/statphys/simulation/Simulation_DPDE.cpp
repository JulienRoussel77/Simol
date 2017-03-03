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

  //------------- DPDE --------------------
  
  void simulate(DPDE& dyna, System& syst)
  {
    //--- compute energy at the beginning of the step ---
    double old_energy = 0;
    double new_mech_energy = 0;
    double new_int_energy = 0;
    if (dyna.doProjectionDPDE())
      {
	for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
	  old_energy += syst(iOfParticle).potentialEnergy() + syst(iOfParticle).kineticEnergy() + syst(iOfParticle).internalEnergy() ;
	//cout << "old = " << old_energy << endl;
      }
    
    //-- first do a certain number of steps of Hamiltonian part, prescribed by "MTS frequency" in the input file
    //   !! works only for NBody systems... !!
    //cout << "time step is in fact... " << dyna.timeStep() << endl;
    for (int k = 0; k < dyna.MTSfrequency(); k++)
      {
	for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
	  dyna.verletFirstPart(syst(iOfParticle));
	syst.computeAllForces();
	for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
	  dyna.verletSecondPart(syst(iOfParticle));
      }
    
    //---------- Isolated systems -------------
    if (syst.name() == "Isolated")
      {
	//-- fluctuation/dissipation --
	for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
	  //dyna.energyReinjection(syst(iOfParticle));  // integration of p at fixed gamma + energy reinjection
	  dyna.metropolizedEnergyReinjection(syst(iOfParticle));  // Metropolis correction of effective dynamics on p 
      }
    //-- then integrate once the fluctuation terms, with effective timestep MTSfrequency()*timeStep_ --
    //------------- NBody --------------------
    else if (syst.name() == "NBody")
      {
	//-- fluctuation/dissipation --
	if (dyna.gamma() > 0)
	  dyna.fluctuationDissipation(syst);
	//-- thermal conduction --
	if (dyna.kappa() > 0)
	  dyna.thermalConduction(syst);
	
      }
    else 
      throw std::runtime_error("The system is not implemented for DPDE dynamics!");
    
  
    //--- adjust energy at the end of the step by rescaling the internal energies ---
    if (dyna.doProjectionDPDE())
      {
	for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
	  {
	    new_mech_energy += syst(iOfParticle).potentialEnergy() + syst(iOfParticle).kineticEnergy();
	    new_int_energy += syst(iOfParticle).internalEnergy();
	  }
	//cout << "new = " << new_mech_energy + new_int_energy << endl;
	double rescaling_factor = (old_energy-new_mech_energy)/new_int_energy;
	for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
	  syst(iOfParticle).internalEnergy() *= rescaling_factor;
      }
  }

}
  
 
