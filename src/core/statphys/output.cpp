#include "output.hpp"

using std::cout; 
using std::endl; 

namespace simol 
{
  

  Output::Output(std::string const& outputFoldername):outputFoldername_(outputFoldername), outParticles_(outputFoldername+"particles.txt"), outReplica_(outputFoldername+"replica.txt"),verbose_(1)
  {
    assert(outParticles_.is_open());
    assert(outReplica_.is_open());
    std::cout << "Output written in " << outputFoldername_ << std::endl;
  }

  void Output::display(Particle const& particle, double time)
  {
    outParticles_ << time 
		  << " " << particle.position() 
		  << " " << particle.momentum() 
		  << " " << particle.kineticEnergy()
		  << " " << particle.potentialEnergy()
		  << " " << particle.energy()
		  << std::endl;
  }
  
  void Output::finalDisplay(Particle const& particle, double time)
  {
    outReplica_ << time 
		  << " " << particle.position() 
		  << " " << particle.momentum() 
		  << " " << particle.kineticEnergy()
		  << " " << particle.potentialEnergy()
		  << " " << particle.energy()
		  << std::endl;
  }
  
}