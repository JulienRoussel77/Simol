#include "output.hpp"

using std::cout; 
using std::endl; 

namespace simol 
{
  
    Statistics::Statistics(int nbIndices):values_(nbIndices, 0), nbValues_(nbIndices, 0){};
    
    void Statistics::append(int i, double value)
    {
      values_[i] = (nbValues_[i] * values_[i] + value) / (++nbValues_[i]);
    }
    
    double& Statistics::operator()(int i)
    {
      return values_[i];
    }
    
    const double& Statistics::operator()(int i) const
    {
      return values_[i];
    }
  
  Output::Output(Input const& input):
    outputFoldername_(input.outputFoldername()), 
    outParticles_(input.outputFoldername()+"particles.txt"), 
    outReplica_(input.outputFoldername()+"replicas.txt"),
    outCorrelation_(input.outputFoldername()+"correlation.txt"),
    verbose_(1),
    timeStep_(input.timeStep()),
    responseForces_(input.dimension()),
    velocityRef_(input.dimension()),
    timeRef_(0),
    decorrelationNumberOfIterations_(input.decorrelationNumberOfIterations()),
    autocorrelationV(decorrelationNumberOfIterations_),
    autocorrelationF(decorrelationNumberOfIterations_)
  {
    assert(outParticles_.is_open());
    assert(outReplica_.is_open());
    std::cout << "Output written in " << input.outputFoldername() << std::endl;
    velocityRef_(0) = input.initialMomentum();
    cout << "decorrelationNumberOfIterations : " << decorrelationNumberOfIterations_ << endl;
  }
  
  void Output::initialize(dvec const& velocityRef)
  {
    responseForces_.fill(0);
    //integratedAutocorrelationP_ = 0;
    velocityRef_ = velocityRef;
  }
  
  int& Output::verbose()
  {
    return verbose_;
  }
  
  const int& Output::verbose() const
  {
    return verbose_;
  }
  
  const double& Output::timeStep() const
  {
    return timeStep_;
  }
  
  dvec& Output::responseForces() 
  {
    return responseForces_;
  }
  
  const dvec& Output::responseForces() const
  {
    return responseForces_;
  }
  
  const double& Output::responseForces(int const& i) const
  {
    return responseForces_(i);
  }
  
  double& Output::responseForces(int const& i)
  {
    return responseForces_(i);
  }
  
  /*double& Output::integratedAutocorrelationP()
  {
    return integratedAutocorrelationP_; 
  }*/

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
  
  
  void Output::finalDisplay(Particle const& particle, dvec const& externalForce, double time)
  {
    outReplica_ << time 
 		  << " " << externalForce(0)   
		  << " " << particle.position() 
		  << " " << particle.momentum() 
		  << " " << responseForces(0)
		  << std::endl;
  }
  
    double& Output::timeRef()
    {
      return timeRef_;
    }
    const double& Output::timeRef() const
    {
      return timeRef_;
    }
    
    dvec& Output::velocityRef()
    {
      return velocityRef_;
    }
    const dvec& Output::velocityRef() const
    {
      return velocityRef_;
    }
    
    dvec& Output::forceRef()
    {
      return forceRef_;
    }
    const dvec& Output::forceRef() const
    {
      return forceRef_;
    }
    
    size_t & Output::decorrelationNumberOfIterations()
    {
      return decorrelationNumberOfIterations_;
    }
    const size_t& Output::decorrelationNumberOfIterations() const
    {
      return decorrelationNumberOfIterations_;
    }
  
  void Output::finalDisplayAutocorrelations()
  {
    double integralV = 0;
    double integralF = 0;
    for (size_t i=0; i < decorrelationNumberOfIterations_; i++)
      outCorrelation_ << i * timeStep() 
		  << " " << autocorrelationV(i)
		  << " " << (integralV += autocorrelationV(i)*timeStep())
		  << " " << autocorrelationF(i)
		  << " " << (integralF += autocorrelationF(i)*timeStep())
		  << std::endl;
  }
  
  void Output::appendAutocorrelationV(dvec const& velocity, double time)
  {
    autocorrelationV.append((time - timeRef_)/timeStep_, velocityRef_.dot(velocity));
  }
  
  void Output::appendAutocorrelationF(dvec const& force, double time)
  {
    autocorrelationF.append((time - timeRef_)/timeStep_, forceRef_.dot(force));
  }
  
}