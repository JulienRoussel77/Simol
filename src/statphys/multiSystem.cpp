#include "multiSystem.hpp"

using std::cout; 
using std::endl; 

namespace simol {
  
 MultiSystem::MultiSystem(Input const& input):
  dimension_(input.dimension()),
  nbOfReplicas_(input.nbOfReplicas()), 
  systemReplicas(nbOfReplicas_),  
  dynamicsReplicas(nbOfReplicas_), 
  output(input),
  //controlVariates_(nbOfReplicas_),
  rng_(input.seed(), input.dimension())
 {
      if (nbOfReplicas_ == 1)
     output.verbose() = 1;
   else
     output.verbose() = 0;
   cout << "verbose = " << output.verbose() << endl;
   

   
   std::cout << "nbOfReplicas : " << nbOfReplicas_ << std::endl;
   for (size_t i=0; i<nbOfReplicas_; i++)
   {
     systemReplicas[i] = createSystem(input, i); 
     dynamicsReplicas[i] = createDynamics(input, i);
     dynamicsReplicas[i]->setRNG(&rng_);
     //controlVariates_[i] = createControlVariate(input, dynamicsReplicas[i]->potential(), i);
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
   extForce(0) = externalForceMin_ + indexReplica * (externalForceMax_ - externalForceMin_) / nbOfReplicas_;
   return extForce;
 }
 
 void MultiSystem::launch(Input const& input)
 {
   for (size_t indexOfReplica=0; indexOfReplica < nbOfReplicas_; indexOfReplica++)
   {
     output.reset(input, &dynamicsReplicas[indexOfReplica]->potential(), dynamicsReplicas[indexOfReplica]->galerkin(), indexOfReplica);
     if ((10*indexOfReplica) % nbOfReplicas_ == 0 && nbOfReplicas_ > 1)
       cout << "-- Simulation " << (100 * indexOfReplica) / nbOfReplicas_ << " % completed ----" << endl;
     //dynamicsReplicas[i]->setExternalForce(externalForce(i));
     //if (i>0)
     //  systemReplicas[i]->particle(0) = systemReplicas[i-1]->particle(0);
     systemReplicas[indexOfReplica]->launch(dynamicsReplicas[indexOfReplica], output);
   }
      
 }
  
}