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
    
  public:
    Simulation(Input& input);
    void launch();
    
  protected:
    int dimension_;
    std::shared_ptr<RNG> rng_;
    S system_;
    D dynamics_;
    Output output_;
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

   //-- when using control variates --
   output_.setControlVariates(input, system_.potential(), dynamics_.galerkin());
  }

  template<class D, class S>
  void Simulation<D,S>::launch()
  {
    launchSimu(dynamics_, system_, output_);
  }

  // --------------- Declaration and implementation of external functions -------------------

  //-- fundamental functions --
  template <class D, class S>
  void sampleSystem(D& dyna, S& syst);
  void samplePositions(Dynamics& dyna, System& syst);
  void sampleMomenta(Dynamics& dyna, System& syst);
  void sampleMomenta(LangevinBase& dyna, System& syst);
  template<class D, class S>
  void launchSimu(D& dyna, S& syst, Output& output);

  void simulate(Dynamics& dyna, System& syst);
  
  template <class D>
  void computeOutput(Dynamics const& dyna, System const& syst, Output& output, int iOfIteration);
  void writeOutput(System const& syst, Output& output, int iOfIteration);
  void writeFinalOutput(Dynamics const& dyna, System const& syst, Output& output);

  //-------------------- specializations for various systems -------------------------
  
  //-- specializations for Isolated --
  template <class D>
  void computeOutput(D const& dyna, Isolated const& syst, Output& output, int iOfIteration);
  template <>
  void computeOutput(const DPDE& dyna, Isolated const& syst, Output& output, int iOfIteration);
  void writeFinalOutput(Hamiltonian const& dyna, Isolated const& syst, Output& output);

  //-- specializations for Chain --
  void sampleMomenta(BoundaryLangevin& dyna, Chain& syst);
  void samplePositions(BoundaryLangevin& dyna, Chain& syst);
  void writeOutput(BoundaryLangevin const& dyna, Chain const& syst, Output& output, int iOfIteration);
  template <class D>
  void computeOutput(D const& dyna, Chain const& syst, Output& output, int iOfIteration);
  template <>
  void computeOutput(BoundaryLangevin const& dyna, Chain const& syst, Output& output, int iOfIteration);
  void simulate(BoundaryLangevin& dyna, Chain& syst);
  // Bichain
  void samplePositions(BoundaryLangevin& dyna, BiChain& syst);
  void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output);
  //TriChain
  void samplePositions(BoundaryLangevin& dyna, TriChain& syst);
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output);
  //-- specializations for NBody --
  void samplePositions(Dynamics& dyna, NBody& syst);
  void simulate(Hamiltonian& dyna, NBody& syst);
  template <class D>
  void computeOutput(D const& /*dyna*/, NBody const& syst, Output& output, int /*iOfIteration*/);
  template <>
  void computeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, int /*iOfIteration*/);
  void writeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, int iOfIteration);
  void writeFinalOutput(Hamiltonian const& dyna, NBody const& syst, Output& output);
  

  // ------------------------ specializations for various dynamics ---------------------------------
  
  //-- DPDE --
  void simulate(DPDE& dyna, System& syst);
  void samplePositions(DPDE& dyna, Isolated& syst);
  void writeOutput(DPDE const& dyna, System const& syst, Output& output, int iOfIteration);

  
  //--------------- control variates -----------------
  
  void updateAllControlVariates(Dynamics const& dyna, System const& syst, Output& output, int iOfIteration);
 
  // Hamiltonian
  void updateAllControlVariates(const Hamiltonian& dyna, System const& syst, Output& output, int iOfIteration);
  void writeOutput(Hamiltonian const& /*dyna*/, System const& syst, Output& output, int iOfIteration);
  
  // Langevin
  void updateAllControlVariates(const Langevin& dyna, System const& syst, Output& output, int iOfIteration);
  void writeOutput(Langevin const& /*dyna*/, System const& syst, Output& output, int iOfIteration);
  
  // Overdamped
  Vector<double> generatorOn(const Overdamped& dyna, const System& syst, const ControlVariate& controlVariate);
  
  //Chain
  void updateAllControlVariates(const BoundaryLangevin& dyna, System const& syst, Output& output, int iOfIteration);


  // -------------------- Template implementation --------------------
  
  template <class D, class S>
  void sampleSystem(D& dyna, S& syst)
  {
    cout << "Initialization of the system...";cout.flush();
    
    sampleMomenta(dyna, syst);
    samplePositions(dyna, syst);

    cout << "Done ! / Thermalization...";cout.flush();

    for (int iOfIteration  =0; iOfIteration < dyna.nbOfThermalIterations(); ++iOfIteration)
    {
      syst.thermalize(dyna);
    }

    cout << "Done ! / BurnIn...";cout.flush();

    for (int iOfIteration  =0; iOfIteration < dyna.nbOfBurnInIterations(); ++iOfIteration)
    {
      simulate(dyna, syst);
    }
    cout << "Done !" << endl;
    
    syst.computeAllForces();
  }
  
  template <class D>
  void computeOutput(D const& dyna, System const& syst, Output& output, int iOfIteration)
  {
    output.kineticEnergy() = 0;
    output.potentialEnergy() = 0;
    //-- compute temperature and kinetic energy --
    for (const auto& particle : syst.configuration())
    {
      output.kineticEnergy() += particle.kineticEnergy();
      output.potentialEnergy() += particle.potentialEnergy();
    }
    // In the case of the trichain we add the potential of the wall interaction
    output.potentialEnergy() += syst.boundaryPotEnergy();
    syst.computeProfile(output, dyna, iOfIteration);
    // use control variates
    updateAllControlVariates(dyna, syst, output, iOfIteration);
  }
  
  // Isolated
  
  template <class D>
  void computeOutput(D const& dyna, Isolated const& syst, Output& output, int iOfIteration)
  {
    //-- compute temperature and kinetic energy --
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    updateAllControlVariates(dyna, syst, output, iOfIteration);
  }
  
  // Chain
  template <class D>
  void computeOutput(D const& dyna, Chain const& syst, Output& output, int iOfIteration)
  {throw std::invalid_argument("computeOutput(D const& dyna, Chain const& syst, Output& output, int iOfIteration) not defined !");}
  
  
  // Nbody
  template <class D>
  void computeOutput(D const& dyna, NBody const& syst, Output& output, int iOfIteration)
  {throw std::invalid_argument("computeOutput(D const& dyna, NBody const& syst, Output& output, int iOfIteration) not defined !");}


  // ------------------------------- MAIN Function ----------------------
  
  template<class D, class S>
  void launchSimu(D& dyna, S& syst, Output& output)
  {
    //---- initialization (including burn-in) -----
    sampleSystem(dyna, syst);

    //---- actual iterations -----
    for (int iOfIteration  =0; iOfIteration < dyna.nbOfIterations(); ++iOfIteration)
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
