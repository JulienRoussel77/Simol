#ifndef SIMOL_SIMULATION_NBODY_HPP
#define SIMOL_SIMULATION_NBODY_HPP

#include "Simulation.hpp"


namespace simol {

  void sampleSystem(Dynamics& dyna, NBody& syst);
  void simulate(Hamiltonian& dyna, NBody& syst);
  void computeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, size_t /*iOfIteration*/);
  void writeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, size_t iOfIteration);
  void writeFinalOutput(Dynamics const& dyna, NBody const& syst, Output& output);
  
};
