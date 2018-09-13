#ifndef SIMOL_SIMULATION_HPP
#define SIMOL_SIMULATION_HPP

#include "simol/Tools.hpp"
#include "simol/system/System.hpp"
#include "simol/system/Isolated.hpp"
#include "simol/system/NBody.hpp"
#include "simol/system/Colloid.hpp"
#include "simol/system/Chain.hpp"
#include "simol/system/Bicolor.hpp"

#include "simol/dynamics/BoundaryLangevin.hpp"
#include "simol/dynamics/Hamiltonian.hpp"
#include "simol/dynamics/Overdamped.hpp"
#include "simol/dynamics/Langevin.hpp"
#include "simol/dynamics/DPDE.hpp"

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
    output_(input)
  {
    
    if (input.dynamicsName() != dynamics_.dynamicsName()) 
      throw std::runtime_error("Dynamics generated incompatible with the input file !");
    system_->rng() = rng_;
    dynamics_.rng() = rng_;

    //-- when using control variates --
    //output_.setControlVariates(input, system_->potential(), dynamics_.galerkin());
  }

  
  
  
  

}

#endif
