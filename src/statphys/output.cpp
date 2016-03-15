#include "output.hpp"

using std::cout; 
using std::endl; 
using std::vector;
using std::sin;
using std::cos;

#include <cmath>

namespace simol{
  
  
  Output::Output(Input const& input):
    outputFolderName_(input.outputFolderName()), 
    outObservables_(input.outputFolderName()+"observables.txt"), 
    outParticles_(input.outputFolderName()+"particles.txt"), 
    outReplica_(input.simuTypeName()+"replicas.txt", std::ofstream::app),
    outFinalFlow_(input.simuTypeName()+"finalFlow.txt", std::ofstream::app),
    outFinalVelocity_(input.simuTypeName()+"finalVelocity.txt", std::ofstream::app),
    outCorrelation_(input.outputFolderName()+"correlations.txt"),
    outVelocities_(input.outputFolderName()+"velocitiesChain.txt"),
    outBeam_(input.outputFolderName()+"beamChain.txt"),
    outVelocitiesCV_(input.outputFolderName()+"velocities.txt"),
    outVelocitiesGenerator_(input.outputFolderName()+"velocitiesGenerator.txt"),
    outForcesCV_(input.outputFolderName()+"forces.txt"),
    outLengthsCV_(input.outputFolderName()+"lengths.txt"),
    outMidFlowCV_(input.outputFolderName()+"midFlow.txt"),
    outMidFlowPT_(input.outputFolderName()+"midFlowPost.txt"),
    outSumFlowCV_(input.outputFolderName()+"sumFlow.txt"),
    outSumFlowPT_(input.outputFolderName()+"sumFlowPost.txt"),
    outProfile_(input.outputFolderName()+"profile.txt"),
    outFinalProfile_(input.outputFolderName()+"finalProfile.txt"),
    //profilePath_(input.outputFolderName()+"profile.txt"),
    verbose_(1),
    periodNbOfIterations_(input.outputPeriodNbOfIterations()),
    profilePeriodNbOfIterations_(input.outputProfilePeriodNbOfIterations()),
    timeStep_(0),
    dimension_(input.dimension()),
    nbOfParticles_(input.nbOfParticles()),
    nbOfIterations_(input.nbOfIterations()),
    decorrelationNbOfIterations_(0),
    nbOfAutocoPts_(input.nbOfAutocoPts()),
    doFinalFlow_(input.doFinalFlow()),
    doFinalVelocity_(input.doFinalVelocity())
  {
    std::cout << "Output written in " << input.outputFolderName() << std::endl;
		std::cout << "Final output written in " << input.simuTypeName() << std::endl;
    assert(outParticles_.is_open());
    assert(outReplica_.is_open());
    
    outObservables_ << "# time kineticEnergy potentialEnergy energy temperature" << endl;
    outParticles_ << "# time index position momentum kineticEnergy potentialEnergy energy force" << endl;
    outReplica_ << "# time dt externalForce position momentum <responseForces> <responseForces2>" << endl;
    outVelocities_ << "# time i=0 i=N/4 i=N/2 i=3N/4 i=N-1" << endl;
    outVelocitiesCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
    outForcesCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
    outLengthsCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
    outMidFlowCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl; 
		outSumFlowCV_ << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl; 
		
		std::ofstream outInput(input.outputFolderName()+"inputFile.txt");
		std::ifstream inInput(input.inputPath());
		outInput << input.inputFlux().rdbuf();
	}
  
  void Output::reset(Input const& input, Potential* potential, Galerkin* galerkin, size_t iOfReplica)
  {
		cout << "nbOfParticles : " << nbOfParticles_ << endl;
    timeStep_ = input.timeStep(iOfReplica);
		nbOfIterations_ = input.nbOfIterations(iOfReplica);
    assert(input.outputPeriodTime(iOfReplica) == 0 || input.outputPeriodTime(iOfReplica) >= timeStep());
    decorrelationNbOfIterations_ = input.decorrelationNbOfIterations(iOfReplica);
    cout << "decorrelationNbOfIterations : " << decorrelationNbOfIterations_ << endl;
    verbose_ = (periodNbOfIterations() > 0);
    velocityCV_ = createControlVariate(input, potential, galerkin, iOfReplica);
    forceCV_ = createControlVariate(input, potential, galerkin, iOfReplica);
    lengthCV_ = createControlVariate(input, potential, galerkin, iOfReplica);
    midFlowCV_ = createControlVariate(input, potential, galerkin, iOfReplica);
		sumFlowCV_ = createControlVariate(input, potential, galerkin, iOfReplica);

		cout << "Initializing the profiles : arrays " << decorrelationNbOfIterations() << " X " << nbOfParticles_ << endl;
		kinTempProfile_ = AutocorrelationStats<double>(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_);
		potTempTopProfile_ = AutocorrelationStats<double>(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_);
		potTempBotProfile_ = AutocorrelationStats<double>(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_);
		bendistProfile_ = AutocorrelationStats<double>(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_);
  	flowProfile_ = AutocorrelationStats<double>(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_);
		energyMidFlow() = 0;
		energySumFlow() = 0;
		
		/*ofstream outMap(input.outputFolderName()+"map.txt");
		
		midFlowCV_->displayMap(outMap);
		
		ofstream outGradQMap(input.outputFolderName()+"gradQMap.txt");
		midFlowCV_->displayGradQMap(outGradQMap);
		
		ofstream outGradPMap(input.outputFolderName()+"gradPMap.txt");
		midFlowCV_->displayGradPMap(outGradPMap);*/
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
    return periodNbOfIterations_ * timeStep();
  }
  
  const size_t& Output::periodNbOfIterations() const
  {
    return periodNbOfIterations_;
  }
  
    const size_t& Output::profilePeriodNbOfIterations() const
  {
    return profilePeriodNbOfIterations_;
  }
  
  bool Output::doOutput(size_t iOfIteration) const
  {
		return (periodNbOfIterations() > 0 && iOfIteration % periodNbOfIterations() == 0);
  }
  
  bool Output::doProfileOutput(size_t iOfIteration) const
  {
		return (profilePeriodNbOfIterations() > 0 && iOfIteration % profilePeriodNbOfIterations() == 0);
  }
  
  const size_t& Output::nbOfParticles() const
  {
		return nbOfParticles_;
	}
	
	const size_t& Output::nbOfIterations() const
	{
		return nbOfIterations_;
	}
	
	double Output::finalTime() const
	{
		return nbOfIterations_ * timeStep_;
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
  
  const double& Output::energyMidFlow() const
  {
    return energyMidFlow_;
  }
  
  double& Output::energyMidFlow()
  {
    return energyMidFlow_;
  }
  
  const double& Output::energySumFlow() const
  {
    return energySumFlow_;
  }
  
  double& Output::energySumFlow()
  {
    return energySumFlow_;
  }
  
  double Output::energy() const
  {
    return kineticEnergy_ + potentialEnergy_;
  }
  
  double Output::temperature() const
  {
    return 2 * kineticEnergy_ / (dimension_ * nbOfParticles_); 
  }
  
  bool Output::doComputeCorrelations() const
  {
    //return doComputeCorrelations_;
    return decorrelationNbOfIterations();
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
  
  ControlVariate* Output::midFlowCV()
  {
    return midFlowCV_;
  }
  
  ControlVariate* Output::sumFlowCV()
  {
    return sumFlowCV_;
  }
  
  
  const double& Output::timeStep() const
  {
    return timeStep_;
  }
  
  double& Output::timeStep()
  {
    return timeStep_;
  }
  
    size_t & Output::decorrelationNbOfIterations()
  {
    return decorrelationNbOfIterations_;
  }
  
  const size_t& Output::decorrelationNbOfIterations() const
  {
    return decorrelationNbOfIterations_;
  }
  
  double Output::decorrelationTime() const
  {
    return decorrelationNbOfIterations() * timeStep();
  }
  
  const int& Output::nbOfAutocoPts() const
  {
		return nbOfAutocoPts_;
	}
	
  double Output::autocoPtsPeriod() const
  {
		return decorrelationTime() / nbOfAutocoPts();
	}

  void Output::display(vector<Particle> const& configuration, size_t iOfIteration)
  {
		//cout << "output : n = " << iOfIteration << endl;
		    
    outObservables_ << iOfIteration * timeStep() 
		    << " " << kineticEnergy()
		    << " " << potentialEnergy()
		    << " " << energy()
		    << " " << temperature()
		    << std::endl;
		
    if (nbOfParticles_ > 1)
      outVelocities_ << iOfIteration * timeStep() 
		     << " " << configuration[0].momentum()
		     << " " << configuration[nbOfParticles_/4].momentum()
		     << " " << configuration[nbOfParticles_/2].momentum()
		     << " " << configuration[3*nbOfParticles_/4].momentum()
		     << " " << configuration[nbOfParticles_-1].momentum()
		     << endl;
		  
    if (nbOfParticles_ > 1)
      outBeam_ << iOfIteration * timeStep() 
		     << " " << configuration[0].position() - 2*configuration[1].position() + configuration[2].position()
				 //<< " " << configuration[0].position() - 2*configuration[1].position() + configuration[2].position()
		     << " " << configuration[(nbOfParticles_-2)/4].position() - 2*configuration[(nbOfParticles_-2)/4+1].position() + configuration[(nbOfParticles_-2)/4+2].position()
		     << " " << configuration[(nbOfParticles_-2)/2].position() - 2*configuration[(nbOfParticles_-2)/2+1].position() + configuration[(nbOfParticles_-2)/2+2].position()
		     << " " << configuration[3*(nbOfParticles_-2)/4].position() - 2*configuration[3*(nbOfParticles_-2)/4+1].position() + configuration[3*(nbOfParticles_-2)/4+2].position()
		     << " " << configuration[nbOfParticles_-3].position() - 2*configuration[nbOfParticles_-2].position() + configuration[nbOfParticles_-1].position()
		     << endl;	   
				 
		     
    velocityCV_->display(outVelocitiesCV_, iOfIteration * timeStep());
    forceCV_->display(outForcesCV_, iOfIteration * timeStep() );
    lengthCV_->display(outLengthsCV_, iOfIteration * timeStep() );
    midFlowCV_->display(outMidFlowCV_, iOfIteration * timeStep() );
		sumFlowCV_->display(outSumFlowCV_, iOfIteration * timeStep() );
		
		//displayGeneratorOnBasis(outVelocitiesGenerator_, configuration, velocityCV_, iOfIteration * timeStep());
		
		if (doProfileOutput(iOfIteration))
		{
			for (size_t i = 0; i < nbOfParticles_; i++)
				outParticles_ << iOfIteration * timeStep() 
		    << " " << i
		    << " " << configuration[i].position() 
		    << " " << configuration[i].momentum() 
		    << " " << configuration[i].kineticEnergy()
		    << " " << configuration[i].potentialEnergy()
		    << " " << configuration[i].energy()
		    << " " << configuration[i].force()		    
		    << endl;
				
			if (nbOfParticles_ > 1)
			{
				/*outProfile_.close();
				outProfile_.open(outputFolderName_ + "/profiles/profile" + to_string(iOfIteration * timeStep()) + ".txt");
				writeProfile(iOfIteration);
				outProfile_.open(outputFolderName_ + "/profile.txt", std::ofstream::app);*/
				writeProfile(outProfile_, iOfIteration);
			}
		}
  }
  
  void Output::displayGeneratorOnBasis(ofstream& out, vector<Particle> const& configuration, ControlVariate* controlVariate, double time)
	{
		//out << time << " " << configuration[0].position(0) << " " << configuration[0].momentum(0) << " " << controlVariate->lastGeneratorOnBasis()(0) << endl;
		out << time << " " << modulo(configuration[0].position(0), -M_PI, M_PI) << " " << configuration[0].momentum(0) << " " << controlVariate->lastGeneratorOnBasis()(0) << " " << controlVariate->basisFunction(configuration) << endl;
		
	}
  
  void Output::writeProfile(ofstream& out_, size_t iOfIteration)
	{
		assert(out_.is_open());
		for (size_t iOfParticle = 0; iOfParticle < nbOfParticles_; iOfParticle++)
			out_ << iOfIteration * timeStep() << " " 
					 << iOfParticle << " " 
					<< bendistProfile_.mean(iOfParticle) << " "
					<< bendistProfile_.standardDeviation(iOfParticle) / sqrt(iOfIteration * timeStep()) << " "
					<< flowProfile_.mean(iOfParticle) << " "
					<< flowProfile_.standardDeviation(iOfParticle) / sqrt(iOfIteration * timeStep()) << " "
					<< kinTempProfile_.mean(iOfParticle) << " "
					<< kinTempProfile_.standardDeviation(iOfParticle) / sqrt(iOfIteration * timeStep()) << " "
					<< potTempTopProfile_.mean(iOfParticle) / potTempBotProfile_.mean(iOfParticle) << " "
					<< endl;
	}
  
  
  void Output::finalDisplay(vector<Particle> const& configuration, Vector<double> const& externalForce)
  {
    cout<< "replica : " << finalTime()
				<< " " << timeStep()
				<< " " << externalForce(0)
				<< " " << configuration[0].position() 
				<< " " << configuration[0].momentum() 
				<< " " << forceCV_->meanObservable()
				<< " " << forceCV_->stdDeviationObservable()
				//<< " " << forceCV_->meanBetterObservable()
				//<< " " << forceCV_->stdDeviationBetterObservable()
				<< " " << midFlowCV_->meanObservable()
				<< " " << midFlowCV_->stdDeviationObservable()
				<< " " << sumFlowCV_->meanObservable()
				<< " " << sumFlowCV_->stdDeviationObservable()
				<< std::endl;
		
    assert(outReplica_.is_open());
    
    outReplica_ << finalTime()
								<< " " << timeStep()
								<< " " << externalForce(0)   
								<< " " << configuration[0].position() 
								<< " " << configuration[0].momentum() 
								<< " " << forceCV_->meanObservable()
								<< " " << forceCV_->stdDeviationObservable()
								//<< " " << forceCV_->meanBetterObservable()
								//<< " " << forceCV_->stdDeviationBetterObservable()
								<< " " << midFlowCV_->meanObservable()
								<< " " << midFlowCV_->stdDeviationObservable()
								<< " " << sumFlowCV_->meanObservable()
								<< " " << sumFlowCV_->stdDeviationObservable()
								<< std::endl;
		
		/*for (size_t iOfParticle = 0; iOfParticle < nbOfParticles_; iOfParticle++)
			outProfile_ << iOfParticle << " " 
									<< kinTempProfile_.mean(iOfParticle) << " "
									<< kinTempProfile_.standardDeviation(iOfParticle) / sqrt(finalTime()) << " "
									<< potTempTopProfile_.mean(iOfParticle) << " "
									<< potTempBotProfile_.mean(iOfParticle) << " "
									<< bendingProfile_.mean(iOfParticle) << " "
									<< bendingProfile_.standardDeviation(iOfParticle) / sqrt(finalTime()) << " "
									<< flowProfile_.mean(iOfParticle) << " "
									<< flowProfile_.standardDeviation(iOfParticle) / sqrt(finalTime())
									<< endl;*/
									
		//outProfile_.close();
		//ofstream outFinalProfile_.open(outputFolderName_ + "/finalProfile.txt");
		writeProfile(outFinalProfile_, nbOfIterations());
									
    midFlowCV_->postTreat(outMidFlowPT_, timeStep());
		sumFlowCV_->postTreat(outSumFlowPT_, timeStep());
  }
  
  void Output::displayFinalFlow(double temperature, double delta_temperature, double tau, double xi)
	{
		//cout << "displayFinalFlow(double temperature, double delta_temperature, double tau)";
		if (doFinalFlow_)
		{
			assert(outFinalFlow_.is_open());
			outFinalFlow_ << std::left << setw(10) << finalTime()
			<< " " << setw(5) << timeStep()
			<< " " << setw(6) << nbOfParticles()
			<< " " << setw(4) << temperature
			<< " " << setw(4) << delta_temperature
			<< " " << setw(4) << tau
			<< " " << setw(6) << xi
			<< " " << setw(12) << midFlowCV_->meanObservable()
			<< " " << setw(12) << midFlowCV_->stdDeviationObservable()
			<< " " << setw(12) << sumFlowCV_->meanObservable()
			<< " " << setw(12) << sumFlowCV_->stdDeviationObservable()
			<< std::endl;
		}
	}
	
  void Output::displayFinalVelocity(double temperature, double externalForce, int nbOfFourier, int nbOfHermite)
	{
		//cout << "displayFinalFlow(double temperature, double delta_temperature, double tau)";
		if (doFinalVelocity_)
		{
			assert(outFinalVelocity_.is_open());
			outFinalVelocity_ << std::left << setw(10) << finalTime()
			<< " " << setw(5) << timeStep()
			<< " " << setw(6) << nbOfParticles()
			<< " " << setw(4) << temperature
			<< " " << setw(6) << externalForce
			<< " " << setw(3) << nbOfFourier
			<< " " << setw(3) << nbOfHermite
			<< " " << setw(12) << velocityCV_->meanObservable()
			<< " " << setw(12) << velocityCV_->stdDeviationObservable()
			<< " " << setw(12) << velocityCV_->meanBetterObservable()
			<< " " << setw(12) << velocityCV_->stdDeviationBetterObservable()
			<< std::endl;
		}
	}

 
  
  void Output::finalDisplayAutocorrelations()
  {
    cout << "decorrelationNbOfIterations : " << decorrelationNbOfIterations_ << endl;
    double integralV = 0;
    double integralF = 0;
    double integralQ = 0;
    double midFlowQ = 0;
		double sumFlowQ = 0;
    //VectorXd integralMidFlowB2(midFlowCV_->nbOfFunctions());
		size_t midNb = nbOfParticles_/2;
    for (int i=0; i < nbOfAutocoPts(); i++)
    {
      outCorrelation_ << i * autocoPtsPeriod() 
		  << " " << velocityCV_->autocorrelation(i) - pow(velocityCV_->meanObservable(), 2)
		  << " " << (integralV += velocityCV_->autocorrelation(i)*autocoPtsPeriod())
		  << " " << forceCV_->autocorrelation(i) - pow(forceCV_->meanObservable(), 2)
		  << " " << (integralF += forceCV_->autocorrelation(i)*autocoPtsPeriod())
		  << " " << lengthCV_->autocorrelation(i) - pow(lengthCV_->meanObservable(), 2)
		  << " " << (integralQ += lengthCV_->autocorrelation(i)*autocoPtsPeriod())
		  << " " << midFlowCV_->autocorrelation(i) - pow(midFlowCV_->meanObservable(), 2)
		  << " " << (midFlowQ += midFlowCV_->autocorrelation(i)*autocoPtsPeriod())
			<< " " << sumFlowCV_->autocorrelation(i) - pow(sumFlowCV_->meanObservable(), 2)
		  << " " << (sumFlowQ += sumFlowCV_->autocorrelation(i)*autocoPtsPeriod());
      //for (size_t iOfFunction = 0; iOfFunction < flowCV_->nbOfFunctions(); iOfFunction++)
		//		outCorrelation_  << " " << flowCV_->autocorrelationB2(i, iOfFunction) - pow(flowCV_->correlationB2(iOfFunction), 2)
		//			<< " " << (integralFlowB2(iOfFunction) += flowCV_->autocorrelationB2(i, iOfFunction)*timeStep());
			outCorrelation_ << " " << bendistProfile_(i, midNb) - pow(bendistProfile_.mean(midNb), 2)
      << std::endl;
    }
  }
 
 	void Output::appendKinTempProfile(double value, size_t iOfIteration, size_t iOfParticle)
	{
		kinTempProfile_.append(value, iOfIteration, iOfParticle);
	}
	
	void Output::appendPotTempTopProfile(double value, size_t iOfIteration, size_t iOfParticle)
	{
		potTempTopProfile_.append(value, iOfIteration, iOfParticle);
	}
	
		 	void Output::appendPotTempBotProfile(double value, size_t iOfIteration, size_t iOfParticle)
	{
		potTempBotProfile_.append(value, iOfIteration, iOfParticle);
	}
	
	void Output::appendBendistProfile(double value, size_t iOfIteration, size_t iOfParticle)
	{
		bendistProfile_.append(value, iOfIteration, iOfParticle);
	}
	
	void Output::appendFlowProfile(double value, size_t iOfIteration, size_t iOfParticle)
	{
		flowProfile_.append(value, iOfIteration, iOfParticle);
	}
  
}
