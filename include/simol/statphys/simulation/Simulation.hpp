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
      // pseudo-random number generator (Mersenne-Twister)
      std::shared_ptr<RNG> rng_;
      // /!\ System is stored as a System* so in all the library: printName(*system_) will return "System" but system_->printName() will return "DaughterSystem"
      System* system_;
      // Dynamics of the system (Hamiltonian, Langevin, Overdamped, ...)
      D dynamics_;
      shared_ptr<CVBasis> cvBasis_;
      Galerkin* galerkin_;
      // object managing the processing of the data and the ouput of relevent information
      Output output_;
      
  };



  // --------------- Declaration and implementation of external functions -------------------
  
  System* createSystem(Input const& input);
  
  
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
