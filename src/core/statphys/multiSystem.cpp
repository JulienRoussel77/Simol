#include "multiSystem.hpp"

using std::cout; 
using std::endl; 

namespace simol {
  
 MultiSystem::MultiSystem(Input const& input):
  size_(input.numberOfReplica()), 
  replica(size_),   
  output(input.outputFoldername())
 {
   std::cout << "numberOfReplica : " << size_ << std::endl;
   for (size_t i=0; i<size_; i++)
     replica[i] = createSystem(input);
   
   if (size_ == 1)
     output.verbose_ = 1;
   else
     output.verbose_ = 0;
   cout << "verbose = " << output.verbose_ << endl;
 }
 
 MultiSystem::~MultiSystem()
 {
   for (auto&& system : replica)
     delete system;
 }
 
 void MultiSystem::launch(Dynamics* model, double const& timeStep, int const& numberOfIterations)
 {
   for (size_t i=0; i<size_; i++)
   {
     if ((10*i) % size_ == 0 && size_ > 1)
       cout << "---- " << (100 * i) / size_ << " % completed ----" << endl;
     replica[i]->launch(model, output, timeStep, numberOfIterations);
   }
 }
  
}