#ifndef SIMOL_DYNAMICS_IPP
#define SIMOL_DYNAMICS_IPP

#include "dynamics.hpp"


namespace simol
{

  Dynamics::Dynamics(YAML::Node const& input)
  : potential_(input["Physics"]["Potential"]["Parameter"].as<double>(), 2*M_PI/input["Geometry"]["Length"].as<size_t>())
  {}
  
  /*Dynamics* createDynamics(Potential const& potential)
  {
    return new Hamiltonian(potential);
  }*/
  
  Dynamics* createDynamics(YAML::Node const& input)
  {
    std::string methodname = input["Physics"]["Model"]["Name"].as<std::string>();
    if (methodname == "Hamiltonian")
      return new Hamiltonian(input);
    else if (methodname == "Langevin")
      return new Langevin(input);
    else
      std::cerr << "Method not valid !" << std::endl;
    
    return 0;
  }

  Potential const & Dynamics::potential() const
  { return potential_; }


  
    Hamiltonian::Hamiltonian(YAML::Node const& input):Dynamics(input)
  {}
  

    Langevin::Langevin(YAML::Node const& input):Dynamics(input)
  {}
}

#endif
