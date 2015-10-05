#include "output.hpp"

namespace simol 
{
  

  Output::Output(std::string const& outputFilename):outputFilename_(outputFilename), out_(outputFilename){}

  void Output::display(double const& time, Particle const& particle)
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