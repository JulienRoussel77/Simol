#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{

  //template <>
  void sampleSystem(DPDE& dyna, NBody& syst)
  {
    cout << " Initialization of the system (NBody, DPDE)..." << endl;
    //--- counters for rejections in Metropolis ---
    dyna.rejectionCount() = 0; 
    dyna.negativeEnergiesCount() = 0;
    dyna.totalCountForRejection() = 0;
    //--- configuration ---
    if (syst.restart())
    {
      cout << " - reading from restart file " << syst.restartFileName() << endl;
      ifstream initialConditions(syst.restartFileName());
      assert(initialConditions.is_open());
      int Dim = syst.dimension();
      for (int i = 0; i < syst.nbOfParticles(); i++)
      {
        for (int dim = 0; dim < Dim; dim++)
          initialConditions >> syst.getParticle(i).position(dim) >> std::ws; 
        for (int dim = 0; dim < Dim; dim++)
          initialConditions >> syst.getParticle(i).momentum(dim) >> std::ws;
        initialConditions >> syst.getParticle(i).internalEnergy() >> std::ws;
      }
    }
    else
    {
      sampleMomenta(dyna, syst);
      samplePositions(dyna, syst);
      sampleInternalEnergies(dyna, syst);
      //--- thermalization ---
      cout << " - Thermalization (" << dyna.thermalizationNbOfSteps() << " steps)..." << endl;
      for (long int iOfStep  = 0; iOfStep < dyna.thermalizationNbOfSteps(); ++iOfStep)
        thermalize(dyna,syst);
    }
    //--- calcul des forces ---
    syst.computeAllForces();
    //--- burn-in ---
    cout << " - Burn-in (" << dyna.burninNbOfSteps() << " steps)..." << endl;
    for (long int iOfStep  = 0; iOfStep < dyna.burninNbOfSteps(); ++iOfStep)
      simulate(dyna, syst);
    // note that at the end of that, some non zero average rejection rate is already obtained
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

  //--- DPDE dynamics ---
  void thermalize(DPDE& dyna, NBody& syst)
  {
    //-- Verlet part --
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.secondPartThermalization(syst(iOfParticle));
  }

  ///
  /// A SUPPRIMER ? NE DEVRAIT PAS DEPENDRE DU SYSTEME...
  void simulate(DPDE& dyna, NBody& syst)
  {
    //-- Verlet part --
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletSecondPart(syst(iOfParticle));
    //-- fluctuation/dissipation --
    syst.fluctuationDissipationDPDE(dyna);
  }


  /*void writeOutput(Hamiltonian const& dyna, NBody const& syst, Output& output, long int iOfStep)
  {
    if (output.doOutput(iOfStep))
      output.displayThermoVariables(iOfStep);
    if (output.doLongPeriodOutput(iOfStep))
      output.displayXMakeMol(syst.configuration(), iOfStep, syst.latticeParameter()*syst.nbOfParticlesPerDimension());
  }
  
  void writeOutput(Overdamped const& dyna, NBody const& syst, Output& output, long int iOfStep)
  {    
    if (output.doOutput(iOfStep))
      output.displayThermoVariables(iOfStep);

    if (output.doLongPeriodOutput(iOfStep))
    {
      if (output.doOutParticles()) output.displayParticles(syst.configuration(), iOfStep);
      if (output.doOutXMakeMol()) output.displayXMakeMol(syst.configuration(), iOfStep);   
      if (output.doOutBackUp()) output.displayBackUp(syst.configuration(), iOfStep);   
    }
  }

  void writeOutput(Langevin const& dyna, NBody const& syst, Output& output, long int iOfStep)
  {
    if (output.doOutput(iOfStep))
      output.displayThermoVariables(iOfStep);
    if (output.doLongPeriodOutput(iOfStep))
      output.displayXMakeMol(syst.configuration(), iOfStep, syst.latticeParameter()*syst.nbOfParticlesPerDimension());
  }

  void writeOutput(DPDE const& dyna, NBody const& syst, Output& output, long int iOfStep)
  {
    if (output.doOutput(iOfStep))
      output.displayThermoVariablesDPDE(syst.configuration(),iOfStep);
    if (output.doLongPeriodOutput(iOfStep))
    {
      output.displayXMakeMol(syst.configuration(), iOfStep, syst.latticeParameter()*syst.nbOfParticlesPerDimension());
      output.displayBackUp(syst.configuration(), iOfStep); 
    }
  }*/

  
}
