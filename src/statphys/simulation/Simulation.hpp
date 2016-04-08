#ifndef SIMOL_SIMULATION_HPP
#define SIMOL_SIMULATION_HPP

#include "Tools.hpp"
#include "System.hpp"
#include "Isolated.hpp"
#include "NBody.hpp"
#include "chain/Chain.hpp"
#include "ControlVariate.hpp"

#include "dynamics/BoundaryLangevin.hpp"
#include "dynamics/Hamiltonian.hpp"
#include "dynamics/Overdamped.hpp"
#include "dynamics/Langevin.hpp"
#include "dynamics/DPDE.hpp"

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

  // Global
  
  void sampleSystem(Dynamics& dyna, System& syst);
  Vector<double> generatorOn(const Dynamics& dyna, const System& syst, const ControlVariate& controlVariate);
  void updateAllControlVariates(const Dynamics& dyna, const System& syst, Output& output, size_t iOfIteration);
  void simulate(Dynamics& dyna, System& syst);
  
  template <class D>
  void computeOutput(Dynamics const& dyna, System const& syst, Output& output, size_t iOfIteration);
  void writeOutput(System const& syst, Output& output, size_t iOfIteration);
  void writeFinalOutput(Dynamics const& dyna, System const& syst, Output& output);
  void updateAllControlVariates(Dynamics const& dyna, System const& syst, Output& output, size_t iOfIteration);


  //Isolated
  void sampleSystem(Hamiltonian& /*dyna*/, Isolated& /*syst*/);
  void sampleSystem(LangevinBase& dyna, Isolated& syst);
  void sampleSystem(Overdamped& dyna, Isolated& syst);
  template <class D>
  void computeOutput(D const& dyna, Isolated const& syst, Output& output, size_t iOfIteration);
  template <>
  void computeOutput(const DPDE& dyna, Isolated const& syst, Output& output, size_t iOfIteration);
  void writeFinalOutput(Hamiltonian const& dyna, Isolated const& syst, Output& output);

  

  //Chains
  void sampleSystem(BoundaryLangevin& dyna, BiChain& syst);
  void sampleSystem(BoundaryLangevin& dyna, TriChain& syst);

  void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output);
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output);

  void simulate(BoundaryLangevin& dyna, Chain& syst);
  void updateAllControlVariates(const BoundaryLangevin& dyna, System const& syst, Output& output, size_t iOfIteration);
  template <class D>
  void computeOutput(D const& dyna, Chain const& syst, Output& output, size_t iOfIteration);
  template <>
  void computeOutput(BoundaryLangevin const& dyna, Chain const& syst, Output& output, size_t iOfIteration);
  void writeOutput(BoundaryLangevin const& dyna, Chain const& syst, Output& output, size_t iOfIteration);
  
  //NBody
  void sampleSystem(Dynamics& dyna, NBody& syst);
  void simulate(Hamiltonian& dyna, NBody& syst);
  template <class D>
  void computeOutput(D const& /*dyna*/, NBody const& syst, Output& output, size_t /*iOfIteration*/);
  template <>
  void computeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, size_t /*iOfIteration*/);
  void writeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, size_t iOfIteration);
  void writeFinalOutput(Hamiltonian const& dyna, NBody const& syst, Output& output);
  
    // ------------------------ Classified by Dynamics ---------------------------------
  

  //Hamiltonian
  //Vector<double> generatorOn(const Hamiltonian& dyna, System const& syst, ControlVariate const& controlVariate);
  void updateAllControlVariates(const Hamiltonian& dyna, System const& syst, Output& output, size_t iOfIteration);
  void writeOutput(Hamiltonian const& /*dyna*/, System const& syst, Output& output, size_t iOfIteration);
  
  //Langevin
  //Vector<double> generatorOn(const Langevin& dyna, System const& syst, const ControlVariate& controlVariate);
  void updateAllControlVariates(const Langevin& dyna, System const& syst, Output& output, size_t iOfIteration);
  void writeOutput(Langevin const& /*dyna*/, System const& syst, Output& output, size_t iOfIteration);
  
  //Overdamped
  Vector<double> generatorOn(const Overdamped& dyna, const System& syst, const ControlVariate& controlVariate);
  
  //DPDE
  void simulate(DPDE& dyna, System& syst);
  void sampleSystem(DPDE& dyna, Isolated& syst);
  
  void writeOutput(DPDE const& dyna, System const& syst, Output& output, size_t iOfIteration);
  
  //General
  template<class D, class S>
  void launchSimu(D& dyna, S& syst, Output& output);









// -------------------- Templates implementation --------------------
  
  // Global
  
  template <class D>
  void computeOutput(D const& dyna, System const& syst, Output& output, size_t iOfIteration)
  {
    output.kineticEnergy() = 0;
    output.potentialEnergy() = 0;
    //Calcul de la température et de l'énergie
    for (const auto& particle : syst.configuration())
    {
      output.kineticEnergy() += particle.kineticEnergy();
      output.potentialEnergy() += particle.potentialEnergy();
    }
    // In the case of the trichain we add the potential of the wall interaction
    output.potentialEnergy() += syst.boundaryPotEnergy();
    syst.computeProfile(output, dyna, iOfIteration);
    updateAllControlVariates(dyna, syst, output, iOfIteration);
  }
  
  // Isolated
  
  template <class D>
  void computeOutput(D const& dyna, Isolated const& syst, Output& output, size_t iOfIteration)
  {
    //cout << "computeOutput(D const& dyna, Isolated const& syst, Output& output, size_t iOfIteration)" << endl;
    //Calcul de la température et de l'énergie
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    updateAllControlVariates(dyna, syst, output, iOfIteration);
  }
  
  // Chain
  
  template <class D>
  void computeOutput(D const& dyna, Chain const& syst, Output& output, size_t iOfIteration)
  {throw std::invalid_argument("computeOutput(D const& dyna, Chain const& syst, Output& output, size_t iOfIteration) not defined !");}
  
  
  // Nbody
  
  template <class D>
  void computeOutput(D const& dyna, NBody const& syst, Output& output, size_t iOfIteration)
  {throw std::invalid_argument("computeOutput(D const& dyna, NBody const& syst, Output& output, size_t iOfIteration) not defined !");}


  // ------------------------------- MAIN Function ----------------------
  
  template<class D, class S>
  void launchSimu(D& dyna, S& syst, Output& output)
  {
    dyna.printName();
    syst.printName();
    cout << "Estimated time : " << 3.5 * syst.nbOfParticles()/1024. * dyna.nbOfIterations() / 1e6 << " hours" << endl;

    //---- initialization (including burn-in) -----
    sampleSystem(dyna, syst);
    syst.computeAllForces();  // TO DO : a mettre dans la fonction d'initialisation...

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
