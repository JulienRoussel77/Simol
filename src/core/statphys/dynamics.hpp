k#ifndef SIMOL_DYNAMICS_HPP
#define SIMOL_DYNAMICS_HPP

#include "potential.hpp"
#include<yaml-cpp/yaml.h>

namespace simol
{
  class Dynamics
  {
    public:
      Dynamics(YAML::Node const& input);

      Potential const & potential() const;
      //friend Dynamics* createDynamics(Potential const& potential);
      friend Dynamics* createDynamics(YAML::Node const& input);

    private:
      //ScalarType mass_;
      Potential potential_;
  };
  
  class Hamiltonian : public Dynamics
  {
  public:
    Hamiltonian(YAML::Node const& input);
  };
  
  
    class Langevin : public Dynamics
  {
  public:
    Langevin(YAML::Node const& input);
  };

}

#include "dynamics.cpp"

#endif
