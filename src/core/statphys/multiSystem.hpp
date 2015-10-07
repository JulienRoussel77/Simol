#pragma once

#include "particleSystem.hpp"

namespace simol {
  
 class MultiSystem
 { 
   size_t size_;
   std::vector<ParticleSystem*> replica;
   Output output;
 public:
   MultiSystem(Input const& input);
   virtual ~MultiSystem();
   void launch(Dynamics* model, double const& timeStep, int const& numberOfIterations);
   int size();
 };
 
}