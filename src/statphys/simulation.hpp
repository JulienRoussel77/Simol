#ifndef SIMOL_SIMULATION_HPP
#define SIMOL_SIMULATION_HPP

#include "tools.hpp"
#include "particleSystem.hpp"
#include "chain/chain.hpp"
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
  
  ParticleSystem* createSystem(Input  const& input);
  Dynamics* createDynamics(Input  const& input);

  
  void initializeMomenta(const Dynamics& dyna, ParticleSystem& syst);
  //void initializeSystem(Dynamics& dyna, ParticleSystem& syst);
  void simulate(Dynamics& dyna, ParticleSystem& syst);
  void updateAllControlVariates(const Dynamics& dyna, const ParticleSystem& syst, Output& output, size_t iOfIteration); 
  Vector<double> generatorOn(const Dynamics& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate);
  void computeOutput(const Dynamics& dyna, const ParticleSystem& syst, Output& output,  size_t iOfIteration);
  //Isolated
  void initializeSystem(Dynamics& dyna, Isolated& syst);
  void initializeSystem(Langevin& dyna, Isolated& syst);
  //Hamiltonian
  Vector<double> generatorOn(const Hamiltonian& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate);
  //Langevin
  Vector<double> generatorOn(const Langevin& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate);
  void updateAllControlVariates(const Langevin& dyna, const ParticleSystem& syst, Output& output, size_t iOfIteration);
  //Overdamped
  Vector<double> generatorOn(const Overdamped& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate);
  //Chains
  void initializeMomenta(const BoundaryLangevin& dyna, Chain& syst);
  void initializeSystem(BoundaryLangevin& dyna, BiChain& syst);
  void initializeSystem(BoundaryLangevin& dyna, TriChain& syst);
  void simulate(BoundaryLangevin& dyna, Chain& syst);
  Vector<double> generatorOn(const BoundaryLangevin& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate); 

  //template<typename D, typename S>
  //void launchSimu(Dynamics& dyna, ParticleSystem& syst, Output& output);
  
  
  
  // Implementation des templates
  

  

}

#endif
