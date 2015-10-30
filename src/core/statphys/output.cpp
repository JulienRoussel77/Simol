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
    outVelocitiesCV_(input.outputFoldername()+"velocities.txt"),
    outForcesCV_(input.outputFoldername()+"forces.txt"),
    outLengthsCV_(input.outputFoldername()+"lengths.txt"),
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
    cout << "decorrelationNumberOfIterations : " << decorrelationNumberOfIterations_ << endl;
    
    outObservables_ << "# time kineticEnergy potentialEnergy energy temperature" << endl;
    outParticles_ << "# time index position momentum kineticEnergy potentialEnergy energy force" << endl;
    outReplica_ << "# time index externalForce position momentum responseForces" << endl;
    outVelocities_ << "# time i=0 i=N/4 i=N/2 i=3N/4 i=N-1" << endl;
    outVelocitiesCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
    outForcesCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
    outLengthsCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
  }
  
  void Output::reset(Input const& input, Potential& potential, size_t indexOfReplica)
  {
    timeStep_ = input.timeStep(indexOfReplica);
    assert(input.outputPeriodTime(indexOfReplica) == 0 || input.outputPeriodTime(indexOfReplica) >= timeStep());
    decorrelationNumberOfIterations_ = input.decorrelationNumberOfIterations(indexOfReplica);
    //autocorrelationVelocity_ = AutocorrelationStats<dvec>(decorrelationNumberOfIterations_, timeStep_);
    //autocorrelationForce_ = AutocorrelationStats<dvec>(decorrelationNumberOfIterations_, timeStep_);
    //length_ = AutocorrelationStats<double>(decorrelationNumberOfIterations_, timeStep_);
    velocityCV_ = createControlVariate(input, potential, indexOfReplica);
    forceCV_ = createControlVariate(input, potential, indexOfReplica);
    lengthCV_ = createControlVariate(input, potential, indexOfReplica);
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
		     
    displayControlVariate(outVelocitiesCV_, velocityCV(), indexOfIteration * timeStep() );
    displayControlVariate(outForcesCV_, forceCV(), indexOfIteration * timeStep() );
    displayControlVariate(outLengthsCV_, lengthCV(), indexOfIteration * timeStep() );
    //displayVelocity(indexOfIteration);
  }
  
  
  void Output::finalDisplay(vector<Particle> const& configuration, dvec const& externalForce, size_t indexOfIteration)
  {
    outReplica_ << " " << timeStep()
		<< " " << timeStep() * indexOfIteration
		<< " " << externalForce(0)   
		<< " " << configuration[0].position() 
		<< " " << configuration[0].momentum() 
		<< " " << forceCV_->meanObservable()
		<< " " << forceCV_->stdDeviationObservable()
		<< std::endl;
  }
  
  /*void Output::displayVelocity(size_t indexOfIteration)
  {
    outVelocities_ << indexOfIteration * timeStep() 
		<< " " << autocorrelationVelocity_.lastValue()
		<< " " << autocorrelationVelocity_.mean()
		<< " " << autocorrelationVelocity_.standardDeviation()
		<< std::endl;
    
  }*/
  
  void Output::displayControlVariate(std::ofstream& out, ControlVariate const* controleVariate, double time) const
  {
    out << time
	<< " " << controleVariate->lastB()
	<< " " << controleVariate->meanB()
	<< " " << controleVariate->stdDeviationB() * sqrt(decorrelationTime())
	<< " " << controleVariate->lastD()
	<< " " << controleVariate->meanD()
	<< " " << controleVariate->stdDeviationD() * sqrt(decorrelationTime())
	<< " " << controleVariate->lastObservable()
	<< " " << controleVariate->meanObservable()
	<< " " << controleVariate->stdDeviationObservable() * sqrt(decorrelationTime())
	<< " " << controleVariate->lastBetterObservable()
	<< " " << controleVariate->meanBetterObservable()
	<< " " << controleVariate->stdDeviationBetterObservable() * sqrt(decorrelationTime())
	<< " " << controleVariate->lastGeneratorOnBasis()
	<< " " << controleVariate->meanGeneratorOnBasis()
	<< " " << controleVariate->stdDeviationGeneratorOnBasis() * sqrt(decorrelationTime())
	<< endl;
  }
  
  void Output::finalDisplayAutocorrelations()
  {
    cout << "decorrelationNumberOfIterations : " << decorrelationNumberOfIterations_ << endl;
    double integralV = 0;
    double integralF = 0;
    double integralQ = 0;
    for (size_t i=0; i < decorrelationNumberOfIterations_; i++)
      outCorrelation_ << i * timeStep() 
		  << " " << velocityCV_->autocorrelation(i)
		  << " " << (integralV += velocityCV_->autocorrelation(i)*timeStep())
		  << " " << forceCV_->autocorrelation(i)
		  << " " << (integralF += forceCV_->autocorrelation(i)*timeStep())
		  << " " << lengthCV_->autocorrelation(i)
		  << " " << (integralQ += lengthCV_->autocorrelation(i)*timeStep())		  
		  << std::endl;
    
  }
 
  
}