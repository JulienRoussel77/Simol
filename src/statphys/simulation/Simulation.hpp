#ifndef SIMOL_SIMULATION_HPP
#define SIMOL_SIMULATION_HPP

#include "Tools.hpp"
#include "System.hpp"
#include "chain/Chain.hpp"
#include "ControlVariate.hpp"

#include "dynamics/BoundaryLangevin.hpp"
#include "dynamics/Hamiltonian.hpp"
#include "dynamics/Overdamped.hpp"
#include "dynamics/Langevin.hpp"

namespace simol {

  template<typename D, typename S>
  class Simulation
  {
    int dimension_;
    std::shared_ptr<RNG> rng_;
    S system_;
    D dynamics_;
    Output output_;
  public:
    Simulation(Input& input);

    void launch();
  };

  template<class D, class S>
  Simulation<D,S>::Simulation(Input& input):
    dimension_(input.dimension()),
    rng_(std::make_shared<RNG>(RNG(input.seed(), input.dimension()))),
    system_(input),
    dynamics_(input),
    output_(input)
 {

   system_.rng() = rng_;
   dynamics_.rng() = rng_;

   output_.setControlVariates(input, system_.potential(), dynamics_.galerkin());
 }

  template<class D, class S>
  void Simulation<D,S>::launch()
  {
    //cout << "Simulation::launch" << endl;
    launchSimu(dynamics_, system_, output_);
  }



  // --------------- Declaration and implementation of external functions -------------------


  //template<class D, class S>
  
  
  //void sampleSystem(Dynamics& /*dyna*/, System& /*syst*/){};
  void simulate(Dynamics& dyna, System& syst);
  void updateAllControlVariates(const Dynamics& dyna, const System& syst, Output& output, size_t iOfIteration);
  Vector<double> generatorOn(const Dynamics& dyna, const System& syst, const ControlVariate& controlVariate);
  void computeOutput(Dynamics const& dyna, System const& syst, Output& output, size_t iOfIteration);
  void writeOutput(System const& syst, Output& output, size_t iOfIteration);
  void writeFinalOutput(Dynamics const& dyna, System const& syst, Output& output);
  void updateAllControlVariates(Dynamics const& dyna, System const& syst, Output& output, size_t iOfIteration);


  //Isolated
  //template <class D>
  //void sampleSystem(Dynamics& dyna, Isolated& syst);
  void sampleSystem(Hamiltonian& /*dyna*/, Isolated& /*syst*/);
  void sampleSystem(UniformStochasticDynamics& dyna, Isolated& syst);
  void sampleSystem(Overdamped& dyna, Isolated& syst);
  void computeOutput(Dynamics const& dyna, const Isolated& syst, Output& output, size_t iOfIteration);
  void writeOutput(System const& syst, Output& output, size_t iOfIteration);

  //void writeFinalOutput(Dynamics const& dyna, Isolated const& syst, Output& output);
  void writeFinalOutput(Hamiltonian const& dyna, Isolated const& syst, Output& output);

  //Hamiltonian
  Vector<double> generatorOn(const Hamiltonian& dyna, System const& syst, ControlVariate const& controlVariate);
  void updateAllControlVariates(const Hamiltonian& dyna, System const& syst, Output& output, size_t iOfIteration);
  void writeOutput(Hamiltonian const& /*dyna*/, System const& syst, Output& output, size_t iOfIteration);
  
  //Langevin
  Vector<double> generatorOn(const Langevin& dyna, System const& syst, const ControlVariate& controlVariate);
  void updateAllControlVariates(const Langevin& dyna, System const& syst, Output& output, size_t iOfIteration);
  void writeOutput(Langevin const& /*dyna*/, System const& syst, Output& output, size_t iOfIteration);
  
  //Overdamped
  Vector<double> generatorOn(const Overdamped& dyna, const System& syst, const ControlVariate& controlVariate);

  //Chains
  void sampleSystem(BoundaryLangevin& dyna, BiChain& syst);
  void sampleSystem(BoundaryLangevin& dyna, TriChain& syst);
  void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output);
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output);

  void simulate(BoundaryLangevin& dyna, Chain& syst);
  void updateAllControlVariates(const BoundaryLangevin& dyna, System const& syst, Output& output, size_t iOfIteration);
  Vector<double> generatorOn(const BoundaryLangevin& dyna, System const& syst, const ControlVariate& controlVariate);
  void writeOutput(BoundaryLangevin const& dyna, System const& syst, Output& output, size_t iOfIteration);
  
  //NBody
  void sampleSystem(Dynamics& dyna, NBody& syst);
  void simulate(Hamiltonian& dyna, NBody& syst);
  void computeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, size_t /*iOfIteration*/);
  void writeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, size_t iOfIteration);
  void writeFinalOutput(Hamiltonian const& dyna, NBody const& syst, Output& output);
  
  //DPDE
  void simulate(DPDE& dyna, System& syst);
  void sampleSystem(DPDE& dyna, Isolated& syst);
  void computeOutput(const DPDE& dyna, const Isolated& syst, Output& output, size_t iOfIteration);
  void writeOutput(DPDE const& dyna, System const& syst, Output& output, size_t iOfIteration);

  //General
  template<class D, class S>
  void launchSimu(D& dyna, S& syst, Output& output);









// -------------------- Templates implementation --------------------

  // ------------------------------- MAIN Function ----------------------
  
  template<class D, class S>
  void launchSimu(D& dyna, S& syst, Output& output)
  {
    dyna.printName();
    syst.printName();
    cout << "Estimated time : " << 3.5 * syst.nbOfParticles()/1024. * dyna.nbOfIterations() / 1e6 << " hours" << endl;

    //---- initialization (including burn-in) -----
    sampleSystem(dyna, syst);
    syst.computeAllForces(dyna);  // TO DO : a mettre dans la fonction d'initialisation...

    //---- actual iterations -----
    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfIterations(); ++iOfIteration)
      {
        //--- display progress every time 10% of simulation elapsed ---
        if ((10*iOfIteration) % dyna.nbOfIterations() == 0)
          cout << "---- Run " << (100 * iOfIteration) / dyna.nbOfIterations() << " % completed ----" << endl;

        //--- write outputs if required ----
        computeOutput(dyna, syst, output, iOfIteration);
        writeOutput(dyna, syst, output, iOfIteration);

        //---- update the system by the numerical integration ---
        simulate(dyna, syst);
      }

    //--- write final outputs ----
    writeFinalOutput(dyna, syst, output);
  }







}

#endif
