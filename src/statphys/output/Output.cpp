#include "Output.hpp"

using std::cout; 
using std::endl; 
using std::vector;
using std::sin;
using std::cos;

// #include <cmath>

namespace simol{
  
  
  Output::Output(Input const& input):
    outputFolderName_(input.outputFolderName()), 
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
    
    cout << endl;
    cout << "####################################################" << endl;
    cout << "#" << endl;
    cout << "# Simulation of : " << endl;
    cout << "#    System    : " << input.systemName() << endl;
    cout << "#    Dynamics  : " << input.dynamicsName() << endl;
    cout << "#    Potential : " << input.potentialName() << endl;
    cout << "#" << endl;   
    cout << "# Number of particles  : " << nbOfParticles_ << endl;
    if (nbOfIterations_<1e6)
      cout << "# Number of iterations : " << nbOfIterations_ << endl;
    else
      cout << "# Number of iterations : " << nbOfIterations_/1e6 << "e6" << endl;
    cout << "# Time step            : " << timeStep_ << endl;
    cout << "#" << endl;
    cout << "####################################################" << endl << endl;
    
    std::cout << "Output written in " << input.outputFolderName() << std::endl;
    std::cout << "Final output written in " << input.simuTypeName() << std::endl;
    if (!outParticles_->is_open())
      throw std::runtime_error("The output file does not exist ! Please add it manually.");
    
    std::ofstream outInput(input.outputFolderName()+"inputFile.txt");
    outInput << input.inputFlux().rdbuf();
    
    assert(input.outputPeriodTime() == 0 || input.outputPeriodTime() >= timeStep());
    decorrelationNbOfIterations_ = input.decorrelationNbOfIterations();
  }
	
  void Output::setControlVariates(Input& input, Potential& potential, Galerkin* galerkin)
  {
    velocityCV_ = createControlVariate(input, potential, galerkin);
    forceCV_ = createControlVariate(input, potential, galerkin);
    lengthCV_ = createControlVariate(input, potential, galerkin);
    midFlowCV_ = createControlVariate(input, potential, galerkin);
    sumFlowCV_ = createControlVariate(input, potential, galerkin);
  }
  
  ofstream & Output::outObservables()
  {return *outObservables_;}
  ofstream & Output::outParticles()
  {return *outParticles_;}
  
  ofstream & Output::outFinalVelocity()
  {return *outFinalVelocity_;}
  ofstream & Output::outCorrelation()
  {return *outCorrelation_;}
  ofstream & Output::outVelocitiesCV()
  {return *outVelocitiesCV_;}
  ofstream & Output::outVelocitiesGenerator()
  {return *outVelocitiesGenerator_;}
  ofstream & Output::outForcesCV()
  {return *outForcesCV_;}
  ofstream & Output::outLengthsCV()
  {return *outLengthsCV_;}
  
  double Output::period() const
  {return periodNbOfIterations_ * timeStep();}
  
  const int& Output::periodNbOfIterations() const
  {return periodNbOfIterations_;}
  
  const int& Output::profilePeriodNbOfIterations() const
  {return profilePeriodNbOfIterations_;}
  
  bool Output::doOutput(int iOfIteration) const
  {
    return (periodNbOfIterations() > 0 && iOfIteration % periodNbOfIterations() == 0);
  }
  
  bool Output::doProfileOutput(int iOfIteration) const
  {
    return (profilePeriodNbOfIterations() > 0 && iOfIteration % profilePeriodNbOfIterations() == 0);
  }
  
  const int& Output::nbOfParticles() const
  {return nbOfParticles_;}
  
  const int& Output::nbOfIterations() const
  {return nbOfIterations_;}
  
  double Output::finalTime() const
  {return nbOfIterations_ * timeStep_;}
  
  const double& Output::kineticEnergy() const
  {return kineticEnergy_;}
  
  double& Output::kineticEnergy()
  {return kineticEnergy_;}
  
  const double& Output::potentialEnergy() const
  {return potentialEnergy_;}
  
  double& Output::potentialEnergy()
  {return potentialEnergy_;}
  
  const double& Output::internalEnergy() const
  {return internalEnergy_;}
  
  double& Output::internalEnergy()
  {return internalEnergy_;}
  
  const double& Output::totalVirial() const
  {return totalVirial_;}
  
  double& Output::totalVirial()
  {return totalVirial_;}
  
  double Output::energy() const
  {return kineticEnergy_ + potentialEnergy_;}
  
  double Output::temperature() const
  {return 2 * kineticEnergy_ / (dimension_ * nbOfParticles_); }
  
  double Output::pressure() const
  {return (2 * kineticEnergy_ + totalVirial_) / (dimension_ * nbOfParticles_ * pow(latticeParameter_,dimension_));}
  
  bool Output::doComputeCorrelations() const
  {return decorrelationNbOfIterations();}
  
  ControlVariate& Output::velocityCV()
  {return *velocityCV_;}
  
  ControlVariate& Output::forceCV()
  {return *forceCV_;}
  
  ControlVariate& Output::lengthCV()
  {return *lengthCV_;}
  
  const double& Output::timeStep() const
  {return timeStep_;}
  
  double& Output::timeStep()
  {return timeStep_;}
  
  int & Output::decorrelationNbOfIterations()
  {return decorrelationNbOfIterations_;}
  
  const int& Output::decorrelationNbOfIterations() const
  {return decorrelationNbOfIterations_;}
  
  double Output::decorrelationTime() const
  {return decorrelationNbOfIterations() * timeStep();}
  
  const int& Output::nbOfAutocoPts() const
  {return nbOfAutocoPts_;}
  
  double Output::autocoPtsPeriod() const
  {return decorrelationTime() / nbOfAutocoPts();}
  
  //----------------- display functions --------------------------
  
  void Output::displayObservables(int iOfIteration)
  {
    outObservables() << iOfIteration * timeStep() 
		    << " " << kineticEnergy()
		    << " " << potentialEnergy()
		    << " " << energy()
		    << " " << temperature()
		    << " " << pressure()
		     << std::endl;
  }
  
  void Output::displayParticles(vector<Particle> const& configuration, int iOfIteration)
  {
    for (int i = 0; i < nbOfParticles_; i++)
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
  

  
  void Output::displayObservablesDPDE(vector<Particle> const& /*configuration*/, int iOfIteration)
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
    int midNb = nbOfParticles_/2;
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
	//for (int iOfFunction = 0; iOfFunction < flowCV_->nbOfFunctions(); iOfFunction++)
	//		outCorrelation_  << " " << flowCV_->autocorrelationB2(i, iOfFunction) - pow(flowCV_->correlationB2(iOfFunction), 2)
	//			<< " " << (integralFlowB2(iOfFunction) += flowCV_->autocorrelationB2(i, iOfFunction)*timeStep());
	outCorrelation() << " " << bendistProfile_(i, midNb) - pow(bendistProfile_.mean(midNb), 2)
			 << std::endl;
      }
  }
 

}
