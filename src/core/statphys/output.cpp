#include "output.hpp"

using std::cout; 
using std::endl; 

namespace simol 
{
  

  Output::Output(std::string const& outputFilename):outputFilename_(outputFilename), outAtoms_(outputFilename), outReplica_("replica.txt")
  {
    assert(outAtoms_.is_open());
    assert(outReplica_.is_open());
    std::cout << "Output written in " << outputFilename_ << std::endl;
  }

  void Output::display(Particle const& particle, double time)
  {
    outAtoms_ << time 
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