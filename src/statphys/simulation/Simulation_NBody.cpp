#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{
  template <>
  void sampleSystem(DPDE& dyna, NBody& syst)
  {
    cout << " Initialization of the system (NBody, DPDE)..." << endl;
    sampleMomenta(dyna, syst);
    samplePositions(dyna, syst);
    sampleInternalEnergies(dyna, syst);
    dyna.rejectionCount() = 0; // rejection rate for Metropolis
    syst.computeAllForces();
    //--- thermalization ---
    cout << " - Thermalization (" << dyna.thermalizationNbOfSteps() << " steps)..." << endl;
    for (long int iOfStep  = 0; iOfStep < dyna.thermalizationNbOfSteps(); ++iOfStep)
      thermalize(dyna,syst);
    //--- end of initialization ---
    cout << endl;
    cout << " Starting production mode" << endl;
    cout << endl;
  }

  void sampleInternalEnergies(DPDE const& /*dyna*/, NBody& syst)
  {
     cout << " - Sampling internal energies..." << endl;
     for (int i = 0; i < syst.nbOfParticles(); i++)
       syst.getParticle(i).internalEnergy() = 1;
  }

  void samplePositions(Dynamics& dyna, NBody& syst)
  {
    cout << " - Sampling the positions..." << endl;
    int Dim = syst.dimension();
    int NbPartDim = syst.nbOfParticlesPerDimension();
    double latticeSize = syst.latticeParameter();

    for (int i = 0; i < syst.nbOfParticles(); i++)
      syst.getParticle(i).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(i).mass());
    //-- initialization on a cubic lattice --
    if (Dim == 2)
    {
      for (int i = 0; i < NbPartDim; i++)
        for (int j = 0; j < NbPartDim; j++)
        {
          syst.getParticle(i * NbPartDim + j).position(0) = i * latticeSize;
          syst.getParticle(i * NbPartDim + j).position(1) = j * latticeSize;
        }
    }
    else if (Dim == 3)
    {
      int NbPartDim2 = NbPartDim * NbPartDim;
      for (int i = 0; i < NbPartDim; i++)
        for (int j = 0; j < NbPartDim; j++)
          for (int k = 0; k < NbPartDim; k++)
          {
            syst.getParticle(i * NbPartDim2 + j * NbPartDim + k).position(0) = i * latticeSize;
            syst.getParticle(i * NbPartDim2 + j * NbPartDim + k).position(1) = j * latticeSize;
            syst.getParticle(i * NbPartDim2 + j * NbPartDim + k).position(2) = k * latticeSize;
          }
    }
    else
    {
      throw std::invalid_argument("sampleSystem: Bad dimension, should be 2 or 3");
    }
  }

  //--- Hamiltonian dynamics ---
  void simulate(Hamiltonian& dyna, NBody& syst)
  {
    for (auto && particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces();
    for (auto && particle : syst.configuration())
      dyna.verletSecondPart(particle);
  }

  //--- Langevin dynamics ---
  void simulate(Langevin& dyna, NBody& syst)
  {
    for (auto && particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces();
    for (auto && particle : syst.configuration())
      dyna.updateAfter(particle);
  }

  //--- DPDE dynamics ---
  void thermalize(DPDE& dyna, NBody& syst)
  {
    //-- Verlet part --
    for (auto && particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces();
    for (auto && particle : syst.configuration())
      dyna.secondPartThermalization(particle);
  }

  void simulate(DPDE& dyna, NBody& syst)
  {
    //-- Verlet part --
    for (auto && particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces();
    for (auto && particle : syst.configuration())
      dyna.verletSecondPart(particle);
    //-- fluctuation/dissipation --
    syst.fluctuationDissipationDPDE(dyna);
  }

  template <>
  void computeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, long int /*iOfStep*/)
  {
    output.kineticEnergy() = 0;
    output.potentialEnergy() = 0;
    output.totalVirial() = 0;
    //Calcul de la température et de l'énergie
  for (const auto & particle : syst.configuration())
    {
      output.kineticEnergy() += particle.kineticEnergy();
      output.potentialEnergy() += particle.potentialEnergy();
      output.totalVirial() += particle.virial();
    }
  }

  template <>
  void computeOutput(Langevin const& /*dyna*/, NBody const& syst, Output& output, long int /*iOfStep*/)
  {
    output.kineticEnergy() = 0;
    output.potentialEnergy() = 0;
    output.totalVirial() = 0;
    //Calcul de la température et de l'énergie
  for (const auto & particle : syst.configuration())
    {
      output.kineticEnergy() += particle.kineticEnergy();
      output.potentialEnergy() += particle.potentialEnergy();
      output.totalVirial() += particle.virial();
    }
  }

  template <>
  void computeOutput(DPDE const& dyna, NBody const& syst, Output& output, long int iOfStep)
  {
    output.kineticEnergy() = 0;
    output.potentialEnergy() = 0;
    output.totalVirial() = 0;
    output.internalEnergy() = 0; 
    output.internalTemperature() = 0;
    //-- computation of instantaneous energies --
    for (const auto & particle : syst.configuration())
      {
	output.kineticEnergy() += particle.kineticEnergy();
	output.potentialEnergy() += particle.potentialEnergy();
	output.totalVirial() += particle.virial();
	output.internalEnergy() += particle.internalEnergy();
	output.internalTemperature() += 1/dyna.internalTemperature(particle.internalEnergy());
      }
    output.internalTemperature() /= syst.nbOfParticles(); 
    //-- averages of observables --
    output.appendKineticEnergy(output.kineticEnergy(), iOfStep);
    output.appendPotentialEnergy(output.potentialEnergy(), iOfStep);
    output.appendInternalEnergy(output.internalEnergy(), iOfStep);
    output.appendInternalTemperature(output.internalTemperature(), iOfStep);
    //-- rejection rate --
    output.rejectionCount() = dyna.rejectionCount()/dyna.totalCountForRejection();
  }

  void writeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, long int iOfStep)
  {
    if (output.doOutput(iOfStep))
      output.displayObservables(iOfStep);
    if (output.doLongPeriodOutput(iOfStep))
      output.displayParticlesXMakeMol(syst.configuration(), iOfStep, syst.latticeParameter()*syst.nbOfParticlesPerDimension());
  }

  void writeOutput(Langevin const& /*dyna*/, NBody const& syst, Output& output, long int iOfStep)
  {
    if (output.doOutput(iOfStep))
      output.displayObservables(iOfStep);
    if (output.doLongPeriodOutput(iOfStep))
      output.displayParticlesXMakeMol(syst.configuration(), iOfStep, syst.latticeParameter()*syst.nbOfParticlesPerDimension());
  }

  void writeOutput(DPDE const& /*dyna*/, NBody const& syst, Output& output, long int iOfStep)
  {
    if (output.doOutput(iOfStep))
      output.displayObservablesDPDE(syst.configuration(),iOfStep);
    if (output.doLongPeriodOutput(iOfStep))
      output.displayParticlesXMakeMol(syst.configuration(), iOfStep, syst.latticeParameter()*syst.nbOfParticlesPerDimension());
  }

  
}
