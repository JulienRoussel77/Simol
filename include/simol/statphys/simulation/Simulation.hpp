#ifndef SIMOL_SIMULATION_HPP
#define SIMOL_SIMULATION_HPP

#include "simol/statphys/Tools.hpp"
#include "simol/statphys/system/System.hpp"
#include "simol/statphys/system/Isolated.hpp"
#include "simol/statphys/system/NBody.hpp"
#include "simol/statphys/system/Chain.hpp"
#include "simol/statphys/controlVariate/ControlVariate.hpp"

#include "simol/statphys/dynamics/BoundaryLangevin.hpp"
#include "simol/statphys/dynamics/Hamiltonian.hpp"
#include "simol/statphys/dynamics/Overdamped.hpp"
#include "simol/statphys/dynamics/Langevin.hpp"
#include "simol/statphys/dynamics/DPDE.hpp"

namespace simol
{

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
  Simulation<D, S>::Simulation(Input& input):
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
  void Simulation<D, S>::launch()
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

  //void simulate(Dynamics& dyna, System& syst);

  template <class D, class S>
  void computeOutput(D const& dyna, S const& syst, Output& output, long int iOfStep);
  //void writeOutput(System const& syst, Output& output, long int iOfStep);
  void writeFinalOutput(Dynamics const& dyna, System const& syst, Output& output);

  //-------------------- specializations for various systems -------------------------

  //-- specializations for Isolated --
  void samplePositions(Dynamics& dyna, Isolated& syst);
  template <class D>
  void computeOutput(D const& dyna, Isolated const& syst, Output& output, long int iOfStep);
  void writeFinalOutput(Hamiltonian const& dyna, Isolated const& syst, Output& output);
  void writeFinalOutput(Langevin const& dyna, Isolated const& syst, Output& output);
  void sampleMomenta(Overdamped& /*dyna*/, Isolated& /*syst*/);
  template <>
  void computeOutput(const DPDE& dyna, Isolated const& syst, Output& output, long int iOfStep);
  void samplePositions(DPDE& dyna, Isolated& syst);
  void simulate(DPDE& dyna, System& syst);
  
  //-- specializations for Chain --
  void sampleMomenta(BoundaryLangevin& dyna, Chain& syst);
  void samplePositions(BoundaryLangevin& dyna, Chain& syst);
  void writeOutput(BoundaryLangevin const& dyna, Chain const& syst, Output& output, long int iOfStep);
  //template <class D>
  //void computeOutput(D const& dyna, Chain const& syst, Output& output, long int iOfStep);
  template <>
  void computeOutput(Overdamped const& dyna, Isolated const& syst, Output& output, long int iOfStep);
  template <class S>
  void computeOutput(BoundaryLangevin const& dyna, S const& syst, Output& output, long int iOfStep);
  void simulate(BoundaryLangevin& dyna, Chain& syst);
  // Bichain
  void samplePositions(BoundaryLangevin& dyna, BiChain& syst);
  void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output);
  //TriChain
  void samplePositions(BoundaryLangevin& dyna, TriChain& syst);
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output);
  
  //-- specializations for NBody --
  template <>
  void sampleSystem(DPDE& dyna, NBody& syst);
  void samplePositions(Dynamics& dyna, NBody& syst);
  void sampleInternalEnergies(DPDE const& /*dyna*/, NBody& syst);
  void simulate(Hamiltonian& dyna, NBody& syst);
  void simulate(Overdamped& dyna, NBody& syst);
  void simulate(Langevin& dyna, NBody& syst);
  void simulate(DPDE& dyna, NBody& syst);
  void thermalize(DPDE& dyna, NBody& syst);
  template <class D>
  void computeOutput(D const& /*dyna*/, NBody const& syst, Output& output, long int /*iOfStep*/);
  template <>
  void computeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, long int /*iOfStep*/);
  template <>
  void computeOutput(Overdamped const& /*dyna*/, NBody const& syst, Output& output, long int /*iOfStep*/);
  template <>
  void computeOutput(Langevin const& /*dyna*/, NBody const& syst, Output& output, long int /*iOfStep*/);
  template <>
  void computeOutput(DPDE const& dyna, NBody const& syst, Output& output, long int iOfStep);
  void writeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, long int iOfStep);
  void writeOutput(Overdamped const& /*dyna*/, NBody const& syst, Output& output, long int iOfStep);
  void writeOutput(Langevin const& /*dyna*/, NBody const& syst, Output& output, long int iOfStep);
  void writeOutput(DPDE const& /*dyna*/, NBody const& syst, Output& output, long int iOfStep);
  
  // ------------------------ specializations for various dynamics ---------------------------------
  
  //--- Overdamped ---
  void simulate(Overdamped& dyna, System& syst);
  
  //--- LangevinBase ---
  void simulate(LangevinBase& dyna, System& syst);
  
  //--- Langevin ---
  void simulate(Langevin& dyna, System& syst);

  //-- DPDE --
  void simulate(DPDE& dyna, System& syst);
  void writeOutput(DPDE const& dyna, System const& syst, Output& output, long int iOfStep);
  void writeFinalOutput(DPDE const& dyna, System const& syst, Output& output);
  


  //--------------- control variates -----------------

  void updateAllControlVariates(Dynamics const& dyna, System const& syst, Output& output, long int iOfStep);

  // Hamiltonian
  void updateAllControlVariates(Hamiltonian const& dyna, System const& syst, Output& output, long int iOfStep);
  void writeOutput(Hamiltonian const& /*dyna*/, System const& syst, Output& output, long int iOfStep);

  // Overdamped
  void updateAllControlVariates(Overdamped const& dyna, System const& syst, Output& output, long int iOfStep);
  void writeOutput(Overdamped const& /*dyna*/, System const& syst, Output& output, long int iOfStep);
  
  // Langevin
  void updateAllControlVariates(Langevin const& dyna, System const& syst, Output& output, long int iOfStep);
  void writeOutput(Langevin const& /*dyna*/, System const& syst, Output& output, long int iOfStep);

  // Overdamped
  Vector<double> generatorOn(Overdamped const& dyna, const System& syst, const ControlVariate& controlVariate);

  // Chain
  void updateAllControlVariates(const BoundaryLangevin& dyna, System const& syst, Output& output, long int iOfStep);


  // -------------------- Template implementation --------------------

  template <class D, class S>
  void sampleSystem(D& dyna, S& syst)
  {
    cout << " Initialization of the system..." << endl;

    sampleMomenta(dyna, syst);
    samplePositions(dyna, syst);

    cout << " - Thermalization (" << dyna.thermalizationNbOfSteps() << " steps)..." << endl;

    for (long int iOfStep  = 0; iOfStep < dyna.thermalizationNbOfSteps(); ++iOfStep)
    {
      syst.thermalize(dyna);
    }

    cout << " - Burn-in (" << dyna.burninNbOfSteps() << " steps)..." << endl;

    for (long int iOfStep  = 0; iOfStep < dyna.burninNbOfSteps(); ++iOfStep)
    {
      simulate(dyna, syst);
    }
    cout << " Starting production mode" << endl;
    cout << endl;

    syst.computeAllForces();
  }

  template <class D, class S>
  void computeOutput(D const& dyna, S const& syst, Output& output, long int iOfStep)
  {
    output.kineticEnergy() = 0;
    output.potentialEnergy() = 0;
    //-- compute temperature and kinetic energy --
    for (const auto & particle : syst.configuration())
    {
      output.kineticEnergy() += particle.kineticEnergy();
      output.potentialEnergy() += particle.potentialEnergy();
    }
    // In the case of the trichain we add the potential of the wall interaction
    output.potentialEnergy() += syst.boundaryPotEnergy();
    syst.computeProfile(output, dyna, iOfStep);
    // use control variates
    updateAllControlVariates(dyna, syst, output, iOfStep);
  }

  // Isolated

  template <class D>
  void computeOutput(D const& dyna, Isolated const& syst, Output& output, long int iOfStep)
  {
    //-- compute temperature and kinetic energy --
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    updateAllControlVariates(dyna, syst, output, iOfStep);
  }

  
  /*template <>
  void computeOutput(Overdamped const& dyna, S const& syst, Output& output, long int iOfStep)
  {
    output.potentialEnergy() = 0;
    //Calcul de la température et de l'énergie
    for (const auto & particle : syst.configuration())
    {
      output.potentialEnergy() += particle.potentialEnergy();
    }
    // In the case of the trichain we add the potential of the wall interaction
    output.potentialEnergy() += syst.boundaryPotEnergy();
    updateAllControlVariates(dyna, syst, output, iOfStep);
  }*/

  // Chain
  template <class S>
  void computeOutput(BoundaryLangevin const& dyna, S const& syst, Output& output, long int iOfStep)
  {
    output.kineticEnergy() = 0;
    output.potentialEnergy() = 0;
    //Calcul de la température et de l'énergie
    for (const auto & particle : syst.configuration())
    {
      output.kineticEnergy() += particle.kineticEnergy();
      output.potentialEnergy() += particle.potentialEnergy();
    }
    // In the case of the trichain we add the potential of the wall interaction
    output.potentialEnergy() += syst.boundaryPotEnergy();
    syst.computeProfile(output, dyna, iOfStep);
    updateAllControlVariates(dyna, syst, output, iOfStep);
  }

  template <class D>
  void computeOutput(D const& dyna, Chain const& syst, Output& output, long int iOfStep)
  {throw std::invalid_argument("computeOutput(D const& dyna, Chain const& syst, Output& output, long int iOfStep) not defined !");}


  // Nbody
  template <class D>
  void computeOutput(D const& dyna, NBody const& syst, Output& output, long int iOfStep)
  {throw std::invalid_argument("computeOutput(D const& dyna, NBody const& syst, Output& output, long int iOfStep) not defined !");}


  // ------------------------------- MAIN Function ----------------------

  template<class D, class S>
  void launchSimu(D& dyna, S& syst, Output& output)
  {
    //---- initialization (including burn-in) -----
    sampleSystem(dyna, syst);

    //---- actual steps -----
    for (long int iOfStep  = 0; iOfStep < dyna.nbOfSteps(); ++iOfStep)
    {
      //--- display progress every time 10% of simulation elapsed ---
      if ((10 * iOfStep) % dyna.nbOfSteps() == 0)
        cout << "---- Run " << (100 * iOfStep) / dyna.nbOfSteps() << " % completed ----" << endl;

      //--- write outputs if required ----
      computeOutput(dyna, syst, output, iOfStep);
      writeOutput(dyna, syst, output, iOfStep);

      //---- update the system by the numerical integration ---
      simulate(dyna, syst);
    }

    //--- write final outputs ----
    writeFinalOutput(dyna, syst, output);
  }


}

#endif
