#include "output.hpp"

namespace simol 
{

  Output::Output(std::string const& outputFilename):outputFilename_(outputFilename){}

  void Output::display(double const& time, Particle const& particle)
  {
    /*std::cout << time 
		  << " " << particle.position() 
		  << " " << particle.momentum() 
		  << " " << particle.kineticEnergy()
		  << " " << particle.potentialEnergy(model->potential())
		  << " " << particle.energy(model->potential())
		  << std::endl;*/
  }	
  
}