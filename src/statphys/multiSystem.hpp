
#ifndef SIMOL_MULTISYSTEM_HPP
#define SIMOL_MULTISYSTEM_HPP

#include "tools.hpp"
#include "particleSystem.hpp"
#include "controlVariate.hpp"

namespace simol {

 class MultiSystem
 {
   int dimension_;
   size_t numberOfReplicas_;
   std::vector<ParticleSystem*> systemReplicas;
   std::vector<Dynamics*> dynamicsReplicas;
   Output output;
   std::vector<ControlVariate*> controlVariates_;
   double externalForceMin_, externalForceMax_;
   RNG rng_;
 public:
   MultiSystem(Input const& input);
   virtual ~MultiSystem();
   Vector<double> externalForce(int const& indexOfReplica);
   void launch(Input const& input);
   int size();
 };

}

#endif
