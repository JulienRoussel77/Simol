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

  ///
  /// Main function
  template<class D, class S>
  void Simulation<D, S>::launch()
  {    
    //---- initialization (including burn-in) -----
    sampleSystem(dynamics_, system_);

    //---- actual steps -----
    for (long int iOfStep  = 0; iOfStep < dynamics_.nbOfSteps(); ++iOfStep)
    {
      //--- display progress every time 10% of simulation elapsed ---
      if ((10 * iOfStep) % dynamics_.nbOfSteps() == 0)
        cout << "---- Run " << (100 * iOfStep) / dynamics_.nbOfSteps() << " % completed ----" << endl;

      //--- write outputs if required ----
      computeOutput(dynamics_, system_, output_, iOfStep);
      writeOutput(dynamics_, system_, output_, iOfStep);

      //---- update the system_em by the numerical integration ---
      simulate(dynamics_, system_);
    }

    //--- write final outputs ----
    writeFinalOutput(dynamics_, system_, output_);
  }

  // --------------- Declaration and implementation of external functions -------------------

  //-- fundamental functions --
  template <class D, class S>
  void sampleSystem(Dynamics& dyna, System& syst);
  void sampleSystem(DPDE& dyna, NBody& syst);

  void samplePositions(Dynamics& dyna, System& syst);
  void samplePositions(Dynamics& dyna, Isolated& syst);
  void samplePositions(DPDE& dyna, Isolated& syst);
  void samplePositions(BoundaryLangevin& dyna, BiChain& syst);
  void samplePositions(BoundaryLangevin& dyna, TriChain& syst);
  void samplePositions(Dynamics& dyna, NBody& syst);

  void sampleMomenta(Dynamics& dyna, System& syst);
  void sampleMomenta(Overdamped& dyna, System& syst);
  void sampleMomenta(LangevinBase& dyna, System& syst);
  void sampleMomenta(BoundaryLangevin& dyna, System& syst);
  
  void sampleInternalEnergies(DPDE const& dyna, NBody& syst);
  
  void thermalize(Dynamics& dyna, System& syst);
  void thermalize(Dynamics& dyna, Isolated& syst);
  void thermalize(Dynamics& model, Chain& syst);
  void thermalize(LangevinBase& dyna, Chain& syst);
  void thermalize(DPDE& dyna, NBody& syst);

  void simulate(DPDE& dyna            , NBody& syst);    // A supprimer ?
  void simulate(Dynamics& dyna        , System& syst);
  void simulate(Overdamped& dyna      , System& syst);
  void simulate(LangevinBase& dyna    , System& syst);
  void simulate(Langevin& dyna        , System& syst);
  void simulate(BoundaryLangevin& dyna, System& syst);
  void simulate(DPDE& dyna            , System& syst);
  

  void computeOutput(Dynamics const& dyna, System const& syst, Output& output, long int iOfStep);
  void computeControlVariate(Dynamics const& dyna, vector<Particle*> const& configuration, Output& output);
  
  void writeOutput(Dynamics const& dyna, System const& syst , Output& output, long int iOfStep);

  void writeFinalOutput(Dynamics const& dyna, System const& syst, Output& output);
  //void writeFinalOutput(Hamiltonian const& dyna, Isolated const& syst, Output& output);
  //void writeFinalOutput(Langevin const& dyna, Isolated const& syst, Output& output);
  //void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output);
  //void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output);
  
  
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
      thermalize(dyna, syst);
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


}

#endif
