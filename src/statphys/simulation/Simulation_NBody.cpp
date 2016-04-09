#include "Simulation.hpp"

namespace simol {
 
  void sampleSystem(Dynamics& dyna, NBody& syst)
  {
    int Dim = syst.dimension();   // PAS SUPER, MAIS SINON PBM DE TYPE POUR COMPARAISON ?
    int NbPartDim = syst.nbOfParticlesPerDimension(); 
    double latticeSize = syst.latticeParameter();
    //-- initialization of the momenta according to a Gaussian distribution --
    for (int i = 0; i < syst.nbOfParticles(); i++)
      syst.getParticle(i).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(i).mass());   
    //-- initialization on a cubic lattice --
    if (Dim == 2) 
      {
      for (int i = 0; i < NbPartDim; i++)
        for (int j = 0; j < NbPartDim; j++)
        {
          syst.getParticle(i*NbPartDim+j).position(0) = i*latticeSize;
          syst.getParticle(i*NbPartDim+j).position(1) = j*latticeSize;
        }
    }
    else if (Dim == 3) 
    {
      int NbPartDim2 = NbPartDim*NbPartDim;
      for (int i = 0; i < NbPartDim; i++)
        for (int j = 0; j < NbPartDim; j++)
          for (int k = 0; k < NbPartDim; k++)
          {
            syst.getParticle(i*NbPartDim2+j*NbPartDim+k).position(0) = i*latticeSize;
            syst.getParticle(i*NbPartDim2+j*NbPartDim+k).position(1) = j*latticeSize;
            syst.getParticle(i*NbPartDim2+j*NbPartDim+k).position(2) = k*latticeSize;
          }
    }
    else 
    {
      throw std::invalid_argument("sampleSystem: Bad dimension, should be 2 or 3");
    }
    //cout << "    VERIFICATION : " << syst.nbOfParticlesPerDimension() << endl;
  }
  
  void simulate(Hamiltonian& dyna, NBody& syst)
  {
    for (auto&& particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces();
    for (auto&& particle : syst.configuration())
      dyna.verletSecondPart(particle);
  }
  
  template <>
  void computeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, int /*iOfIteration*/)
  {
    output.kineticEnergy() = 0;
    output.potentialEnergy() = 0;
    output.totalVirial() = 0;
    //Calcul de la température et de l'énergie
    for (const auto& particle : syst.configuration())
      {
      output.kineticEnergy() += particle.kineticEnergy();
      output.potentialEnergy() += particle.potentialEnergy();       
      output.totalVirial() += particle.virial();
      }
  }
  
  void writeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, int iOfIteration)
  {
    if (output.doOutput(iOfIteration))
      output.displayObservables(iOfIteration);
    if (output.doProfileOutput(iOfIteration))
      output.displayParticlesXMakeMol(syst.configuration(), iOfIteration, syst.latticeParameter()*syst.nbOfParticlesPerDimension());        
  }
  
  void writeFinalOutput(Hamiltonian const& /*dyna*/, NBody const& /*syst*/, Output& /*output*/)
  {
    //output.finalDisplay(syst.configuration(), dyna.externalForce());
  }

}
