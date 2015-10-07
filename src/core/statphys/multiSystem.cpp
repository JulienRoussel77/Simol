#include "multiSystem.hpp"

using std::cout; 
using std::endl; 

namespace simol {
  
 MultiSystem::MultiSystem(Input const& input):
  size_(input.numberOfReplica()), 
  replica(size_),   
  output(input.outputFilename())
 {
   std::cout << "numberOfReplica : " << size_ << std::endl;
   for (size_t i=0; i<size_; i++)
     replica[i] = createSystem(input);
 }
 
 MultiSystem::~MultiSystem()
 {
   for (auto&& system : replica)
     delete system;
 }
 
 void MultiSystem::launch(Dynamics* model, double const& timeStep, int const& numberOfIterations)
 {
   for (auto&& system : replica)
   {
     system->launch(model, output, timeStep, numberOfIterations);
   }
 }
  
}