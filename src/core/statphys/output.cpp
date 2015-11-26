#include "output.hpp"

using std::cout; 
using std::endl; 
using std::vector;
using std::sin;
using std::cos;

#include <cmath>

namespace simol{
  
  
  Output::Output(Input const& input):
    outputFoldername_(input.outputFoldername()), 
    outObservables_(input.outputFoldername()+"observables.txt"), 
    outParticles_(input.outputFoldername()+"particles.txt"), 
    outReplica_(input.outputFoldername()+"replicas.txt", std::ofstream::app),
    outCorrelation_(input.outputFoldername()+"correlations.txt"),
    outVelocities_(input.outputFoldername()+"velocitiesChain.txt"),
    outBeam_(input.outputFoldername()+"beamChain.txt"),
    outVelocitiesCV_(input.outputFoldername()+"velocities.txt"),
    outForcesCV_(input.outputFoldername()+"forces.txt"),
    outLengthsCV_(input.outputFoldername()+"lengths.txt"),
    outFlowCV_(input.outputFoldername()+"flow.txt"),
    outFlowPT_(input.outputFoldername()+"flowPost.txt"),
    verbose_(1),
    periodNumberOfIterations_(input.outputPeriodNumberOfIterations()),
    timeStep_(0),
    dimension_(input.dimension()),
    numberOfParticles_(input.numberOfParticles()),
    decorrelationNumberOfIterations_(0)
  {
    std::cout << "Output written in " << input.outputFoldername() << std::endl;
    assert(outParticles_.is_open());
    assert(outReplica_.is_open());
    
    
    outObservables_ << "# time kineticEnergy potentialEnergy energy temperature" << endl;
    outParticles_ << "# time index position momentum kineticEnergy potentialEnergy energy force" << endl;
    outReplica_ << "# time index externalForce position momentum responseForces" << endl;
    outVelocities_ << "# time i=0 i=N/4 i=N/2 i=3N/4 i=N-1" << endl;
    outVelocitiesCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
    outForcesCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
    outLengthsCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
    outFlowCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl; 
  }
  
  void Output::reset(Input const& input, Potential& potential, size_t indexOfReplica)
  {
    timeStep_ = input.timeStep(indexOfReplica);
    assert(input.outputPeriodTime(indexOfReplica) == 0 || input.outputPeriodTime(indexOfReplica) >= timeStep());
    decorrelationNumberOfIterations_ = input.decorrelationNumberOfIterations(indexOfReplica);
    cout << "decorrelationNumberOfIterations : " << decorrelationNumberOfIterations_ << endl;
    verbose_ = (periodNumberOfIterations() > 0);
    velocityCV_ = createControlVariate(input, potential, indexOfReplica);
    forceCV_ = createControlVariate(input, potential, indexOfReplica);
    lengthCV_ = createControlVariate(input, potential, indexOfReplica);
    flowCV_ = createControlVariate(input, potential, indexOfReplica);
  }
    
  int& Output::verbose()
  {
    return verbose_;
  }
  
  const int& Output::verbose() const
  {
    return verbose_;
  }
  
  double Output::period() const
  {
    return periodNumberOfIterations_ * timeStep();
  }
  
  const size_t& Output::periodNumberOfIterations() const
  {
    return periodNumberOfIterations_;
  }
  
  bool Output::doOutput(size_t indexOfIteration) const
  {
    return (periodNumberOfIterations() > 0 && indexOfIteration % periodNumberOfIterations() == 0);
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
  
  const double& Output::energyFlow() const
  {
    return energyFlow_;
  }
  
  double& Output::energyFlow()
  {
    return energyFlow_;
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
  
  ControlVariate* Output::velocityCV()
  {
    return velocityCV_;
  }
  
  ControlVariate* Output::forceCV()
  {
    return forceCV_;
  }
  
  ControlVariate* Output::lengthCV()
  {
    return lengthCV_;
  }
  
  ControlVariate* Output::flowCV()
  {
    return flowCV_;
  }
  
  
  
  const double& Output::timeStep() const
  {
    return timeStep_;
  }
  
  double& Output::timeStep()
  {
    return timeStep_;
  }
  
    size_t & Output::decorrelationNumberOfIterations()
  {
    return decorrelationNumberOfIterations_;
  }
  
  const size_t& Output::decorrelationNumberOfIterations() const
  {
    return decorrelationNumberOfIterations_;
  }
  
  double Output::decorrelationTime() const
  {
    return decorrelationNumberOfIterations() * timeStep();
  }

  void Output::display(vector<Particle> const& configuration, size_t indexOfIteration)
  {
    for (size_t i = 0; i < numberOfParticles_; i++)
      outParticles_ << indexOfIteration * timeStep() 
		    << " " << i
		    << " " << configuration[i].position() 
		    << " " << configuration[i].momentum() 
		    << " " << configuration[i].kineticEnergy()
		    << " " << configuration[i].potentialEnergy()
		    << " " << configuration[i].energy()
		    << " " << configuration[i].force()		    
		    << std::endl;
		    
    outObservables_ << indexOfIteration * timeStep() 
		    << " " << kineticEnergy()
		    << " " << potentialEnergy()
		    << " " << energy()
		    << " " << temperature()
		    << std::endl;
		
    if (numberOfParticles_ > 1)
      outVelocities_ << indexOfIteration * timeStep() 
		     << " " << configuration[0].momentum()
		     << " " << configuration[numberOfParticles_/4].momentum()
		     << " " << configuration[numberOfParticles_/2].momentum()
		     << " " << configuration[3*numberOfParticles_/4].momentum()
		     << " " << configuration[numberOfParticles_-1].momentum()
		     << endl;
		  
    if (numberOfParticles_ > 1)
      outBeam_ << indexOfIteration * timeStep() 
		     << " " << configuration[0].position() - 2*configuration[1].position() + configuration[2].position()
		     << " " << configuration[(numberOfParticles_-2)/4].position() - 2*configuration[(numberOfParticles_-2)/4+1].position() + configuration[(numberOfParticles_-2)/4+2].position()
		     << " " << configuration[(numberOfParticles_-2)/2].position() - 2*configuration[(numberOfParticles_-2)/2+1].position() + configuration[(numberOfParticles_-2)/2+2].position()
		     << " " << configuration[3*(numberOfParticles_-2)/4].position() - 2*configuration[3*(numberOfParticles_-2)/4+1].position() + configuration[3*(numberOfParticles_-2)/4+2].position()
		     << " " << configuration[numberOfParticles_-3].position() - 2*configuration[numberOfParticles_-2].position() + configuration[numberOfParticles_-1].position()
		     << endl;	     
		     
    velocityCV_->display(outVelocitiesCV_, indexOfIteration * timeStep());
    forceCV_->display(outForcesCV_, indexOfIteration * timeStep() );
    lengthCV_->display(outLengthsCV_, indexOfIteration * timeStep() );
    flowCV_->display(outFlowCV_, indexOfIteration * timeStep() );
  }
  
  
  void Output::finalDisplay(vector<Particle> const& configuration, dvec const& externalForce, size_t indexOfIteration)
  {
    cout<< "replica : " << timeStep() * indexOfIteration
		<< " " << timeStep()
		<< " " << externalForce(0)   
		<< " " << configuration[0].position() 
		<< " " << configuration[0].momentum() 
		<< " " << forceCV_->meanObservable()
		<< " " << forceCV_->stdDeviationObservable()
		<< " " << forceCV_->meanBetterObservable()
		<< " " << forceCV_->stdDeviationBetterObservable()
		<< std::endl;
		
    assert(outReplica_.is_open());
    
    outReplica_ << " " << timeStep() * indexOfIteration
		<< " " << timeStep()
		<< " " << externalForce(0)   
		<< " " << configuration[0].position() 
		<< " " << configuration[0].momentum() 
		<< " " << forceCV_->meanObservable()
		<< " " << forceCV_->stdDeviationObservable()
		<< " " << forceCV_->meanBetterObservable()
		<< " " << forceCV_->stdDeviationBetterObservable()
		<< std::endl;
		
    flowCV_->postTreat(outFlowPT_, timeStep());
  }
  

 
  
  void Output::finalDisplayAutocorrelations()
  {
    cout << "decorrelationNumberOfIterations : " << decorrelationNumberOfIterations_ << endl;
    double integralV = 0;
    double integralF = 0;
    double integralQ = 0;
    double flowQ = 0;
    VectorXd integralFlowB2(flowCV_->nbOfFunctions());
    for (size_t i=0; i < decorrelationNumberOfIterations_; i++)
    {
      outCorrelation_ << i * timeStep() 
		  << " " << velocityCV_->autocorrelation(i) - pow(velocityCV_->meanObservable(), 2)
		  << " " << (integralV += velocityCV_->autocorrelation(i)*timeStep())
		  << " " << forceCV_->autocorrelation(i) - pow(forceCV_->meanObservable(), 2)
		  << " " << (integralF += forceCV_->autocorrelation(i)*timeStep())
		  << " " << lengthCV_->autocorrelation(i) - pow(lengthCV_->meanObservable(), 2)
		  << " " << (integralQ += lengthCV_->autocorrelation(i)*timeStep())
		  << " " << flowCV_->autocorrelation(i) - pow(flowCV_->meanObservable(), 2)
		  << " " << (flowQ += flowCV_->autocorrelation(i)*timeStep());
      for (size_t iOfFunction = 0; iOfFunction < flowCV_->nbOfFunctions(); iOfFunction++)
	outCorrelation_  << " " << flowCV_->autocorrelationB2(i, iOfFunction) - pow(flowCV_->correlationB2(iOfFunction), 2)
		  << " " << (integralFlowB2(iOfFunction) += flowCV_->autocorrelationB2(i, iOfFunction)*timeStep());	
      outCorrelation_ << std::endl;
    }
  }
 
  
}