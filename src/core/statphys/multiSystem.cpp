#include "multiSystem.hpp"

using std::cout; 
using std::endl; 

namespace simol {
  
 MultiSystem::MultiSystem(Input const& input):
  dimension_(input.dimension()),
  numberOfReplicas_(input.numberOfReplicas()), 
  systemReplicas(numberOfReplicas_),  
  dynamicsReplicas(numberOfReplicas_), 
  output(input),
  rng_(1, input.dimension())
 {
      if (numberOfReplicas_ == 1)
     output.verbose() = 1;
   else
     output.verbose() = 0;
   cout << "verbose = " << output.verbose() << endl;
   

   
   std::cout << "numberOfReplicas : " << numberOfReplicas_ << std::endl;
   for (size_t i=0; i<numberOfReplicas_; i++)
   {
     systemReplicas[i] = createSystem(input, i); 
     dynamicsReplicas[i] = createDynamics(input, i);
     dynamicsReplicas[i]->setRNG(&rng_);
   }
 }
 
 MultiSystem::~MultiSystem()
 {
   for (auto&& system : systemReplicas)
     delete system;
   for (auto&& dynamics : dynamicsReplicas)
     delete dynamics;
 }
 
 dvec MultiSystem::externalForce(int const& indexReplica){
   dvec extForce(dimension_);
   extForce(0) = externalForceMin_ + indexReplica * (externalForceMax_ - externalForceMin_) / numberOfReplicas_;
   return extForce;
 }
 
 void MultiSystem::launch()
 {
   for (size_t i=0; i<numberOfReplicas_; i++)
   {
     if ((10*i) % numberOfReplicas_ == 0 && numberOfReplicas_ > 1)
       cout << "---- " << (100 * i) / numberOfReplicas_ << " % completed ----" << endl;
     //dynamicsReplicas[i]->setExternalForce(externalForce(i));
     //if (i>0)
     //  systemReplicas[i]->particle(0) = systemReplicas[i-1]->particle(0);
     systemReplicas[i]->launch(dynamicsReplicas[i], output);
   }
      
 }
  
}