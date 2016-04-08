#include "Output.hpp"

using std::cout; 
using std::endl; 
using std::vector;
using std::sin;
using std::cos;

#include <cmath>

namespace simol{
  
  
  Output::Output(Input const& input):
    outputFolderName_(input.outputFolderName()), 
    /*outObservables_(input.outputFolderName()+"observables.txt"), 
    outParticles_(input.outputFolderName()+"particles.txt"), 
    outParticlesXMakeMol_(input.outputFolderName()+"xmakemol.xyz"), 
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
    outFinalProfile_(input.outputFolderName()+"finalProfile.txt"),*/
    periodNbOfIterations_(input.outputPeriodNbOfIterations()),
    profilePeriodNbOfIterations_(input.outputProfilePeriodNbOfIterations()),
    timeStep_(input.timeStep()),
    dimension_(input.dimension()),
    nbOfParticles_(input.nbOfParticles()),
    nbOfIterations_(input.nbOfIterations()),
    latticeParameter_(input.latticeParameter()),
    kineticEnergy_(0),
    potentialEnergy_(0),
    internalEnergy_(0),
    totalVirial_(0),
    energyMidFlow_(0),
    energySumFlow_(0),
    decorrelationNbOfIterations_(input.decorrelationNbOfIterations()),
    nbOfAutocoPts_(input.nbOfAutocoPts()),
    doFinalFlow_(input.doFinalFlow()),
    doFinalVelocity_(input.doFinalVelocity()),
    velocityCV_(nullptr),
    forceCV_(nullptr),
    lengthCV_(nullptr),
    midFlowCV_(nullptr),
    sumFlowCV_(nullptr),
    kinTempProfile_(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_),
    potTempTopProfile_(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_),
    potTempBotProfile_(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_),
    bendistProfile_(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_),
    flowProfile_(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_)
  {
    outObservables_       = std::make_shared<ofstream>(input.outputFolderName()+"observables.txt");
    outObservables() << "# time kineticEnergy potentialEnergy energy temperature" << endl;
    
    outParticles_         = std::make_shared<ofstream>(input.outputFolderName()+"particles.txt");
    outParticles() << "# time index position momentum kineticEnergy potentialEnergy energy force" << endl;

    outParticlesXMakeMol_ = std::make_shared<ofstream>(input.outputFolderName()+"xmakemol.xyz");
    
    outCorrelation_       = std::make_shared<ofstream>(input.outputFolderName()+"correlations.txt");
    outVelocitiesCV_      = std::make_shared<ofstream>(input.outputFolderName()+"velocities.txt");
    outVelocitiesCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;

    outVelocitiesGenerator_ = std::make_shared<ofstream>(input.outputFolderName()+"velocitiesGenerator.txt");
    outForcesCV_          = std::make_shared<ofstream>(input.outputFolderName()+"forces.txt");
    outForcesCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;

    outLengthsCV_         = std::make_shared<ofstream>(input.outputFolderName()+"lengths.txt");
    outLengthsCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
  
    if (doFinalVelocity_)
      outFinalVelocity_     = std::make_shared<ofstream>(input.simuTypeName()+"finalVelocity.txt", std::ofstream::app);

    if (input.systemName() == "BiChain" || input.systemName() == "TriChain")
    {
      outFinalFlow_         = std::make_shared<ofstream>(input.simuTypeName()+"finalFlow.txt", std::ofstream::app);
      outBeam_              = std::make_shared<ofstream>(input.outputFolderName()+"beamChain.txt");
      
      outChainVelocities_   = std::make_shared<ofstream>(input.outputFolderName()+"velocitiesChain.txt");
      outChainVelocities() << "# time i=0 i=N/4 i=N/2 i=3N/4 i=N-1" << endl;
    
      outMidFlowCV_         = std::make_shared<ofstream>(input.outputFolderName()+"midFlow.txt");
      outMidFlowCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl; 

      outMidFlowPT_         = std::make_shared<ofstream>(input.outputFolderName()+"midFlowPost.txt");
      outSumFlowCV_         = std::make_shared<ofstream>(input.outputFolderName()+"sumFlow.txt");
      outSumFlowCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl; 

      outSumFlowPT_         = std::make_shared<ofstream>(input.outputFolderName()+"sumFlowPost.txt");
      outProfile_           = std::make_shared<ofstream>(input.outputFolderName()+"profile.txt");
      
      if (doFinalFlow_)
        outFinalProfile_      = std::make_shared<ofstream>(input.outputFolderName()+"finalProfile.txt");
    }
    
    std::cout << "Output written in " << input.outputFolderName() << std::endl;
    std::cout << "Final output written in " << input.simuTypeName() << std::endl;
    assert(outParticles_->is_open());
    
    std::ofstream outInput(input.outputFolderName()+"inputFile.txt");
    outInput << input.inputFlux().rdbuf();
    
    cout << "nbOfParticles : " << nbOfParticles_ << endl;
    assert(input.outputPeriodTime() == 0 || input.outputPeriodTime() >= timeStep());
    decorrelationNbOfIterations_ = input.decorrelationNbOfIterations();
    cout << "decorrelationNbOfIterations : " << decorrelationNbOfIterations_ << endl;
    
    cout << "Initializing the profiles : arrays " << decorrelationNbOfIterations() << " X " << nbOfParticles_ << endl;
  }
	
	void Output::setControlVariates(Input& input, Potential& potential, Galerkin* galerkin)
  {
    velocityCV_ = createControlVariate(input, potential, galerkin);
    cout << "cv test : " << velocityCV_->nbOfFunctions() << endl;
    forceCV_ = createControlVariate(input, potential, galerkin);
    lengthCV_ = createControlVariate(input, potential, galerkin);
    midFlowCV_ = createControlVariate(input, potential, galerkin);
    sumFlowCV_ = createControlVariate(input, potential, galerkin);
  }
  
  ofstream & Output::outObservables()
  {return *outObservables_;}
  ofstream & Output::outParticles()
  {return *outParticles_;}
  ofstream & Output::outParticlesXMakeMol()
  {return *outParticlesXMakeMol_;}
  ofstream & Output::outFinalFlow()
  {return *outFinalFlow_;}
  ofstream & Output::outFinalVelocity()
  {return *outFinalVelocity_;}
  ofstream & Output::outCorrelation()
  {return *outCorrelation_;}
  ofstream & Output::outChainVelocities()
  {return *outChainVelocities_;}
  ofstream & Output::outBeam()
  {return *outBeam_;}
  ofstream & Output::outVelocitiesCV()
  {return *outVelocitiesCV_;}
  ofstream & Output::outVelocitiesGenerator()
  {return *outVelocitiesGenerator_;}
  ofstream & Output::outForcesCV()
  {return *outForcesCV_;}
  ofstream & Output::outLengthsCV()
  {return *outLengthsCV_;}
  ofstream & Output::outMidFlowCV()    
  {return *outMidFlowCV_;}
  ofstream & Output::outMidFlowPT()
  {return *outMidFlowPT_;}
  ofstream & Output::outSumFlowCV()    
  {return *outSumFlowCV_;}
  ofstream & Output::outSumFlowPT()
  {return *outSumFlowCV_;}
  ofstream & Output::outProfile()
  {return *outProfile_;}
  ofstream & Output::outFinalProfile()
  {return *outFinalProfile_;}
  
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

  const double& Output::internalEnergy() const
  {
    return internalEnergy_;
  }
  
  double& Output::internalEnergy()
  {
    return internalEnergy_;
  }
  
  const double& Output::totalVirial() const
  {
    return totalVirial_;
  }
  
  double& Output::totalVirial()
  {
    return totalVirial_;
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
  
  double Output::pressure() const
  {
    return (2 * kineticEnergy_ + totalVirial_) / (dimension_ * nbOfParticles_ * pow(latticeParameter_,dimension_)); 
  }
  
  bool Output::doComputeCorrelations() const
  {
    //return doComputeCorrelations_;
    return decorrelationNbOfIterations();
  }
  
  ControlVariate& Output::velocityCV()
  {
    return *velocityCV_;
  }
  
  ControlVariate& Output::forceCV()
  {
    return *forceCV_;
  }
  
  ControlVariate& Output::lengthCV()
  {
    return *lengthCV_;
  }
  
  ControlVariate& Output::midFlowCV()
  {
    return *midFlowCV_;
  }
  
  ControlVariate& Output::sumFlowCV()
  {
    return *sumFlowCV_;
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
  
  void Output::displayObservables(size_t iOfIteration)
  {
    outObservables() << iOfIteration * timeStep() 
		    << " " << kineticEnergy()
		    << " " << potentialEnergy()
		    << " " << energy()
		    << " " << temperature()
		    << " " << pressure()
		    << std::endl;
  }
  
  void Output::displayChainPositions(vector<Particle> const& configuration, size_t iOfIteration)
  {
    outBeam() << iOfIteration * timeStep() 
        << " " << configuration[0].position() - 2*configuration[1].position() + configuration[2].position()
      //<< " " << configuration[0].position() - 2*configuration[1].position() + configuration[2].position()
        << " " << configuration[(nbOfParticles_-2)/4].position() - 2*configuration[(nbOfParticles_-2)/4+1].position() + configuration[(nbOfParticles_-2)/4+2].position()
        << " " << configuration[(nbOfParticles_-2)/2].position() - 2*configuration[(nbOfParticles_-2)/2+1].position() + configuration[(nbOfParticles_-2)/2+2].position()
        << " " << configuration[3*(nbOfParticles_-2)/4].position() - 2*configuration[3*(nbOfParticles_-2)/4+1].position() + configuration[3*(nbOfParticles_-2)/4+2].position()
        << " " << configuration[nbOfParticles_-3].position() - 2*configuration[nbOfParticles_-2].position() + configuration[nbOfParticles_-1].position()
        << endl;   
  }
  
  void Output::displayChainMomenta(vector<Particle> const& configuration, size_t iOfIteration)
  {
    outChainVelocities() << iOfIteration * timeStep()
      << " " << configuration[0].momentum()
      << " " << configuration[nbOfParticles_/4].momentum()
      << " " << configuration[nbOfParticles_/2].momentum()
      << " " << configuration[3*nbOfParticles_/4].momentum()
      << " " << configuration[nbOfParticles_-1].momentum()
      << endl;
  }
  
  void Output::displayParticles(vector<Particle> const& configuration, size_t iOfIteration)
  {
    for (size_t i = 0; i < nbOfParticles_; i++)
      outParticles() << iOfIteration * timeStep() 
		    << " " << i
		    << " " << configuration[i].position() 
		    << " " << configuration[i].momentum() 
		    << " " << configuration[i].kineticEnergy()
		    << " " << configuration[i].potentialEnergy()
		    << " " << configuration[i].energy()
		    << " " << configuration[i].force()        
		    << endl;
  }

  //--- display the current configuration in XMakemol format ; specific for NBody systems ---
  void Output::displayParticlesXMakeMol(vector<Particle> const& configuration, size_t iOfIteration, double domainSize)
  {
    outParticlesXMakeMol() << nbOfParticles_ << endl;
    outParticlesXMakeMol() << "Time = " << iOfIteration * timeStep() << endl;
    double coordinate = 0;
    int Dim = dimension_;
    for (size_t i = 0; i < nbOfParticles_; i++)
    {
      outParticlesXMakeMol() << " O  ";
      if (Dim == 3)
      {
	outParticlesXMakeMol_ << " O  ";
	if (Dim == 3)
	{
	  for (int dim = 0; dim < Dim; dim++)
	    {
	      //-- recenter all the coordinates in the interval [-domainSize/2, domainSize/2] --
	      coordinate = configuration[i].position(dim);
	      coordinate -= rint(coordinate/domainSize)*domainSize;
	      outParticlesXMakeMol_ << coordinate << " "; 
	    }
	}
	else if (Dim == 2)
	{
	  for (int dim = 0; dim < Dim; dim++)
	    {
	      //-- recenter all the coordinates in the interval [-domainSize/2, domainSize/2] --
	      coordinate = configuration[i].position(dim);
	      coordinate -= rint(coordinate/domainSize)*domainSize;
	      outParticlesXMakeMol_ << coordinate << " "; 
	    }
	   outParticlesXMakeMol_ << 0 << " ";  
	}
	outParticlesXMakeMol_ << endl;
        for (int dim = 0; dim < dim; dim++)
          {
            //-- recenter all the coordinates in the interval [-domainSize/2, domainSize/2] --
            coordinate = configuration[i].position(dim);
            coordinate -= rint(coordinate/domainSize)*domainSize;
            outParticlesXMakeMol() << coordinate << " "; 
          }
      }
      else if (Dim == 2)
      {
        for (int dim = 0; dim < Dim; dim++)
          {
            //-- recenter all the coordinates in the interval [-domainSize/2, domainSize/2] --
            coordinate = configuration[i].position(dim);
            coordinate -= rint(coordinate/domainSize)*domainSize;
            outParticlesXMakeMol() << coordinate << " "; 
          }
        outParticlesXMakeMol() << 0 << " ";  
      }
      outParticlesXMakeMol() << endl;
    }
  }
  
  void Output::displayProfile(size_t iOfIteration)
  {
    writeProfile(outProfile(), iOfIteration);
  }
  
  void Output::displayObservablesDPDE(vector<Particle> const& /*configuration*/, size_t iOfIteration)
  {
    // garder le vecteur des particules pour evt sortir les positions, etc
    //cout << "output : n = " << iOfIteration << endl;
    double totalEnergy = kineticEnergy() + potentialEnergy() + internalEnergy();
    outObservables() << iOfIteration * timeStep() 
		    << " " << kineticEnergy()
		    << " " << potentialEnergy()
		    << " " << internalEnergy()
		    << " " << totalEnergy 
		    << std::endl;
  }
  

  void Output::displayGeneratorOnBasis(ofstream& out, vector<Particle> const& configuration, ControlVariate& controlVariate, double time)
	{
    out << time << " " << modulo(configuration[0].position(0), -M_PI, M_PI) << " " << configuration[0].momentum(0) << " " << controlVariate.lastGeneratorOnBasis()(0) << " " << controlVariate.basisFunction(configuration) << endl;
		
	}
  
  void Output::writeProfile(ofstream & out_, size_t iOfIteration)
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
  
  
  void Output::finalDisplay(vector<Particle> const& /*configuration*/, Vector<double> const& /*externalForce*/)
  {    
		cout << "Output::finalDisplay" << endl;
		writeProfile(outFinalProfile(), nbOfIterations());
									
    midFlowCV_->postTreat(outMidFlowPT(), timeStep());
		sumFlowCV_->postTreat(outSumFlowPT(), timeStep());
  }
  
  void Output::displayFinalFlow(double temperature, double delta_temperature, double tau, double xi)
	{
    cout << "outFinalFlow_ : " <<  std::left << setw(10) << finalTime()
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
      
		//cout << "displayFinalFlow(double temperature, double delta_temperature, double tau)";
		if (doFinalFlow_)
		{
			assert(outFinalFlow().is_open());
			outFinalFlow() << std::left << setw(10) << finalTime()
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
			assert(outFinalVelocity().is_open());
			outFinalVelocity() << std::left << setw(10) << finalTime()
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
    double integralV = 0;
    double integralF = 0;
    double integralQ = 0;
    double midFlowQ = 0;
		double sumFlowQ = 0;
		size_t midNb = nbOfParticles_/2;
    for (int i=0; i < nbOfAutocoPts(); i++)
    {
      outCorrelation() << i * autocoPtsPeriod() 
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
			outCorrelation() << " " << bendistProfile_(i, midNb) - pow(bendistProfile_.mean(midNb), 2)
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
