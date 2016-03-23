#ifndef SIMOL_SIMULATION_HPP
#define SIMOL_SIMULATION_HPP

#include "tools.hpp"
#include "particleSystem.hpp"
#include "controlVariate.hpp"

namespace simol {

  class Simulation
  {
    int dimension_;
    ParticleSystem* system_;
    Dynamics* dynamics_;
    Output output_;
    //ControlVariate* controlVariate_;
    RNG rng_;
  public:
    Simulation(Input& input);
    virtual ~Simulation(); 
 
    void launch();
  };

  void initializeMomenta(const Dynamics& dyna, ParticleSystem& syst);
  void initializeSystem(Dynamics& dyna, ParticleSystem& syst);
  void simulate(Dynamics& dyna, ParticleSystem& syst);
  void updateAllControlVariates(const Dynamics& dyna, const ParticleSystem& syst, Output& output, size_t iOfIteration); 
  DenseMatrix<double> generatorOn(const Dynamics& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate);
  void computeOutput(const Dynamics& dyna, const ParticleSystem& syst, Output& output,  size_t iOfIteration);
  //Hamiltonian
  DenseMatrix<double> generatorOn(const Hamiltonian& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate);
  //Langevin
  DenseMatrix<double> generatorOn(const Langevin& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate);
  void updateAllControlVariates(const Langevin& dyna, const ParticleSystem& syst, Output& output, size_t iOfIteration);
  //Overdamped
  DenseMatrix<double> generatorOn(const Overdamped& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate);
  //Chains
  void initializeMomenta(const BoundaryLangevin& dyna, Chain& syst);
  void initializeSystem(BoundaryLangevin& dyna, BiChain& syst);
  void initializeSystem(BoundaryLangevin& dyna, TriChain& syst);
  void simulate(BoundaryLangevin& dyna, Chain& syst);
  DenseMatrix<double> generatorOn(const BoundaryLangevin& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate); 

  void launchSimu(Dynamics& dyna, ParticleSystem& syst, Output& output);
  

}

#endif
