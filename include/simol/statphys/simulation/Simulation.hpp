#ifndef SIMOL_SIMULATION_HPP
#define SIMOL_SIMULATION_HPP

#include "simol/statphys/Tools.hpp"
#include "simol/statphys/system/System.hpp"
#include "simol/statphys/system/Isolated.hpp"
#include "simol/statphys/system/NBody.hpp"
#include "simol/statphys/system/Colloid.hpp"
#include "simol/statphys/system/Chain.hpp"
#include "simol/statphys/system/Bicolor.hpp"
#include "simol/statphys/controlVariate/ControlVariate.hpp"

#include "simol/statphys/dynamics/BoundaryLangevin.hpp"
#include "simol/statphys/dynamics/Hamiltonian.hpp"
#include "simol/statphys/dynamics/Overdamped.hpp"
#include "simol/statphys/dynamics/Langevin.hpp"
#include "simol/statphys/dynamics/DPDE.hpp"

namespace simol
{

  template<typename D>
  class Simulation
  {

    public:
      Simulation(Input& input);
      void launch();

    public:
      int dimension_;
      std::shared_ptr<RNG> rng_;
      // /!\ System is stored as a System* so in all the library: printName(*system_) will return "System" but system_->printName() will return "DaughterSystem"
      System* system_;
      D dynamics_;
      shared_ptr<CVBasis> cvBasis_;
      Galerkin* galerkin_;
      Output output_;
      
  };



  // --------------- Declaration and implementation of external functions -------------------
  
  System* createSystem(Input const& input);

  /*//-- fundamental functions --
  void sampleSystem(Dynamics& dyna, System& syst);
  
  void sampleInternalEnergies(Dynamics const& dyna, System& syst);
  void sampleInternalEnergies(DPDE const& dyna, System& syst);
  
  void thermalize(Dynamics& dyna, System& syst);
  void thermalize(BoundaryLangevin& dyna, System& syst);
  void thermalize(DPDE& dyna, System& syst);
  
  void computeOutput(Dynamics const& dyna, System const& syst, Output& output, long int iOfStep);
  void computeControlVariate(System const& syst, Output& output);
  void writeOutput(Dynamics const& dyna, System const& syst , Output& output, long int iOfStep);
  void writeFinalOutput(System const& syst, Output& output);*/
  
  
  // -------------------- Template implementation --------------------
  
  template<class D>
  Simulation<D>::Simulation(Input& input):
    dimension_(input.dimension()),
    rng_(std::make_shared<RNG>(RNG(input.seed(), input.dimension()))),
    system_(createSystem(input)),
    dynamics_(input),
    cvBasis_(createCVBasis(input)),
    //cvBasis_(input, make_shared<DynamicsParameters>(dynamics_.parameters())),
    galerkin_(createGalerkin(input, cvBasis_)),
    output_(input, cvBasis_)
  {
    if (cvBasis_ && !cvBasis_->cvCoeffs_)
      throw runtime_error("cvBasis_.cvCoeffs_ not initialized !");
    
    if (input.dynamicsName() != dynamics_.dynamicsName()) 
      throw std::runtime_error("Dynamics generated incompatible with the input file !");
    system_->rng() = rng_;
    dynamics_.rng() = rng_;

    //-- when using control variates --
    //output_.setControlVariates(input, system_->potential(), dynamics_.galerkin());
  }

  
  
  
  

}

#endif
