#include "output.hpp"

using std::cout; 
using std::endl; 

namespace simol 
{
  
  Output::Output(Input const& input):
    outputFoldername_(input.outputFoldername()), 
    outParticles_(input.outputFoldername()+"particles.txt"), 
    outReplica_(input.outputFoldername()+"replicas.txt"),
    verbose_(1),
    sumForces_(input.dimension())
  {
    assert(outParticles_.is_open());
    assert(outReplica_.is_open());
    std::cout << "Output written in " << input.outputFoldername() << std::endl;
  }
  
  void Output::initialize()
  {
    sumForces_.fill(0);
  }
  
  dvec& Output::sumForces() 
  {
    return sumForces_;
  }
  
  std::vector<dvec>& Output::responseForces()
  {
    return responseForces_;
  }
  
  const dvec& Output::responseForces(int const& i) const
  {
    return responseForces_[i];
  }
  
  dvec& Output::responseForces(int const& i)
  {
    return responseForces_[i];
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
  
  void Output::finalDisplayExternalForce(Particle const& particle, dvec const& externalForce, dvec const& responseForce, double time)
  {
    outReplica_ << time 
 		  << " " << externalForce(0)   
		  << " " << particle.position() 
		  << " " << particle.momentum() 
		  << " " << responseForce(0)

		  << std::endl;
  }
  
}