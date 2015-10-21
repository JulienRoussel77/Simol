#include "output.hpp"

using std::cout; 
using std::endl; 
using std::vector;

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
    outObservables_(input.outputFoldername()+"observables.txt"), 
    outParticles_(input.outputFoldername()+"particles.txt"), 
    outReplica_(input.outputFoldername()+"replicas.txt"),
    outCorrelation_(input.outputFoldername()+"correlation.txt"),
    verbose_(1),
    timeStep_(input.timeStep()),
    dimension_(input.dimension()),
    numberOfParticles_(input.numberOfParticles()),
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
  
  const double& Output::kineticEnergy() const
  {
    return kineticEnergy_;
  }
  
  double& Output::kineticEnergy()
  {
    return kineticEnergy_;
  }
  
  const double& Output::potentialEnergy() const
  {
    return potentialEnergy_;
  }
  
  double& Output::potentialEnergy()
  {
    return potentialEnergy_;
  }
  
  double Output::energy() const
  {
    return kineticEnergy_ + potentialEnergy_;
  }
  
  double Output::temperature() const
  {
    return 2 * kineticEnergy_ / (dimension_ * numberOfParticles_); 
  }
  
  bool Output::doComputeCorrelations() const
  {
    //return doComputeCorrelations_;
    return decorrelationNumberOfIterations();
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

  void Output::display(vector<Particle> const& configuration, double time)
  {
    for (size_t i = 0; i < configuration.size(); i++)
      outParticles_ << time 
		    << " " << i
		    << " " << configuration[i].position() 
		    << " " << configuration[i].momentum() 
		    << " " << configuration[i].kineticEnergy()
		    << " " << configuration[i].potentialEnergy()
		    << " " << configuration[i].energy()
		    << " " << configuration[i].force()		    
		    << std::endl;
		    
    outObservables_ << time
		    << " " << kineticEnergy()
		    << " " << potentialEnergy()
		    << " " << energy()
		    << " " << temperature()
		    << std::endl;
    
  }
  
  
  void Output::finalDisplay(vector<Particle> const& configuration, dvec const& externalForce, double time)
  {
    for (size_t i = 0; i < configuration.size(); i++)
      outReplica_ << time 
		  << " " << i
 		  << " " << externalForce(0)   
		  << " " << configuration[i].position() 
		  << " " << configuration[i].momentum() 
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