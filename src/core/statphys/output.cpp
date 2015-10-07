#include "output.hpp"

using std::cout; 
using std::endl; 

namespace simol 
{
  

  Output::Output(std::string const& outputFilename):outputFilename_(outputFilename), out_(outputFilename)
  {
    assert(out_.is_open());
    std::cout << "Output written in " << outputFilename_ << std::endl;
  }

  void Output::display(Particle const& particle, double time)
  {
    out_ << time 
		  << " " << particle.position() 
		  << " " << particle.momentum() 
		  << " " << particle.kineticEnergy()
		  << " " << particle.potentialEnergy()
		  << " " << particle.energy()
		  << std::endl;
  }	
  
}