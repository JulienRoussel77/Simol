#include "Output.hpp"

using std::cout; 
using std::endl; 
using std::vector;

namespace simol{
  
  ///
  /// creation of output class
  /// including creation/opening of appropriate output files
  Output::Output(Input const& input):
    outputFolderName_(input.outputFolderName()), 
    printPeriodNbOfSteps_(input.printPeriodNbOfSteps()),
    printLongPeriodNbOfSteps_(input.printLongPeriodNbOfSteps()),
    timeStep_(input.timeStep()),
    dimension_(input.dimension()),
    nbOfParticles_(input.nbOfParticles()),
    nbOfSteps_(input.nbOfSteps()),
    latticeParameter_(input.latticeParameter()),
    kineticEnergy_(0),
    potentialEnergy_(0),
    internalEnergy_(0),
    totalVirial_(0),
    energyMidFlow_(0),
    energySumFlow_(0),
    decorrelationNbOfSteps_(input.decorrelationNbOfSteps()),
    nbOfAutocoPts_(input.nbOfAutocoPts()),
    doFinalFlow_(input.doFinalFlow()),
    doFinalVelocity_(input.doFinalVelocity()),
    velocityCV_(nullptr),
    forceCV_(nullptr),
    lengthCV_(nullptr),
    midFlowCV_(nullptr),
    sumFlowCV_(nullptr),
    kinTempProfile_(decorrelationNbOfSteps(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_),
    potTempTopProfile_(decorrelationNbOfSteps(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_),
    potTempBotProfile_(decorrelationNbOfSteps(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_),
    bendistProfile_(decorrelationNbOfSteps(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_),
    flowProfile_(decorrelationNbOfSteps(), decorrelationTime(), nbOfAutocoPts(), nbOfParticles_)
  {
    //-- standard observables in this file --
    outObservables_       = std::make_shared<ofstream>(input.outputFolderName()+"observables.txt");
    outObservables() << "# time kineticEnergy potentialEnergy energy temperature" << endl;
    
    //-- longer outputs if required, e.g. configuration of the system --
    if (printLongPeriodNbOfSteps_ > 0)
    { 
      if (input.systemName() == "NBody")
      {
        outParticlesXMakeMol_ = std::make_shared<ofstream>(input.outputFolderName()+"xmakemol.xyz");
      }
      else
      {
        outParticles_ = std::make_shared<ofstream>(input.outputFolderName()+"particles.txt");
        outParticles() << "# time index position momentum kineticEnergy potentialEnergy energy force" << endl;
      }
    }

    //-- for autocorrelations --
    if (decorrelationNbOfSteps_ > 0)
      outCorrelation_ = std::make_shared<ofstream>(input.outputFolderName()+"correlations.txt");
    
    //-- average velocities for 1D nonequilibrium --
    if ( (input.systemName() == "Isolated") && (input.dynamicsName() == "Langevin") )
    {
      outVelocitiesCV_      = std::make_shared<ofstream>(input.outputFolderName()+"velocities.txt");
      outVelocitiesCV() << "# time b <b> <b2> D <D> <D2> observable <observable> <observable2> LPhi <LPhi> <LPhi2>" << endl;  
    }
    
    //-- average forces for 1D nonequilibrium overdamped --
    if ( (input.systemName() == "Isolated") && (input.dynamicsName() == "Overdamped") )
    {  
      outForcesCV_ = std::make_shared<ofstream>(input.outputFolderName()+"forces.txt");
      outForcesCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
    }
    
    //-- used for control variates --
    if (input.controlVariateName() != "None")
      outVelocitiesGenerator_ = std::make_shared<ofstream>(input.outputFolderName()+"velocitiesGenerator.txt");
    
    //---- TO DO: check status ---
    outLengthsCV_         = std::make_shared<ofstream>(input.outputFolderName()+"lengths.txt");
    outLengthsCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
  
    //--- to have a summary of the fina, estimated average velocity ---
    if (doFinalVelocity_)
      outFinalVelocity_     = std::make_shared<ofstream>(input.simuTypeName()+"finalVelocity.txt", std::ofstream::app);

    //-- outputs specific to chains --
    if (input.systemName() == "BiChain" || input.systemName() == "TriChain")
    {
      outFinalProfile_      = std::make_shared<ofstream>(input.outputFolderName()+"finalProfile.txt");
      
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
        outFinalFlow_         = std::make_shared<ofstream>(input.simuTypeName()+"finalFlow.txt", std::ofstream::app);
    }
    
    //--------- announce where input/outputs are ----------------------------
    std::cout << endl;
    std::cout << "Output written in " << input.outputFolderName() << std::endl;
    std::cout << "Final output written in " << input.simuTypeName() << std::endl;
    
    if (outParticles_ && !outParticles_->is_open())
      throw std::runtime_error("The output file does not exist ! Please add it manually.");
    
    //-- copy read input into file --
    std::ofstream outInput(input.outputFolderName()+"inputFile.txt");
    outInput << input.inputFlux().rdbuf();
      
    //------------------ screen output for control --------------------------
    cout << endl;
    cout << "-----------------------------------------------------" << endl;
    cout << endl;
    cout << " Simulation of : " << endl;
    cout << "    System    : " << input.systemName() << endl;
    cout << "    Dynamics  : " << input.dynamicsName() << endl;
    cout << "    Potential : " << input.potentialName() << endl;
    cout << endl;   
    cout << " Number of particles  : " << nbOfParticles_ << endl;
    if (nbOfSteps_<1e6)
      cout << " Number of steps : " << nbOfSteps_ << endl;
    else
      cout << " Number of steps : " << nbOfSteps_/1e6 << "e6" << endl;
    cout << " Time step            : " << timeStep_ << endl;
    cout << "" << endl;
    cout << "-----------------------------------------------------" << endl << endl;
    
  }
	
	//----- various assessors --------
	
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
  
  double Output::printPeriodTime() const
  {return printPeriodNbOfSteps_ * timeStep();}
  
  const int& Output::printPeriodNbOfSteps() const
  {return printPeriodNbOfSteps_;}
  
  double Output::printLongPeriodTime() const
  {return printLongPeriodNbOfSteps_ * timeStep();}
  
  const int& Output::printLongPeriodNbOfSteps() const
  {return printLongPeriodNbOfSteps_;}
  
  const int& Output::nbOfParticles() const
  {return nbOfParticles_;}
  
  const int& Output::nbOfSteps() const
  {return nbOfSteps_;}
  
  double Output::finalTime() const
  {return nbOfSteps_ * timeStep_;}
  
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
  {return decorrelationNbOfSteps();}
  
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
  
  int & Output::decorrelationNbOfSteps()
  {return decorrelationNbOfSteps_;}
  
  const int& Output::decorrelationNbOfSteps() const
  {return decorrelationNbOfSteps_;}
  
  double Output::decorrelationTime() const
  {return decorrelationNbOfSteps() * timeStep();}
  
  const int& Output::nbOfAutocoPts() const
  {return nbOfAutocoPts_;}
  
  double Output::autocoPtsPeriod() const
  {return decorrelationTime() / nbOfAutocoPts();}
  
  //---------------- check whether outputs should be performed --------------------
  
  bool Output::doOutput(int iOfStep) const
  {
    return (printPeriodNbOfSteps() > 0 && iOfStep % printPeriodNbOfSteps() == 0);
  }
  
  bool Output::doLongPeriodOutput(int iOfStep) const
  {
    return (printLongPeriodNbOfSteps() > 0 && iOfStep % printLongPeriodNbOfSteps() == 0);
  }
  
  //----------------- display functions --------------------------
  
  void Output::displayObservables(int iOfStep)
  {
    outObservables() << iOfStep * timeStep() 
		    << " " << kineticEnergy()
		    << " " << potentialEnergy()
		    << " " << energy()
		    << " " << temperature()
		    << " " << pressure()
		     << std::endl;
  }
  
  void Output::displayParticles(vector<Particle> const& configuration, int iOfStep)
  {
    for (int i = 0; i < nbOfParticles_; i++)
      outParticles() << iOfStep * timeStep() 
		     << " " << i
		     << " " << configuration[i].position() 
		     << " " << configuration[i].momentum() 
		     << " " << configuration[i].kineticEnergy()
		    << " " << configuration[i].potentialEnergy()
		     << " " << configuration[i].energy()
		     << " " << configuration[i].force()        
		     << endl;
  }

    void Output::displayFinalVelocity(double temperature, double externalForce, int nbOfFourier, int nbOfHermite)
  {
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
  
  void Output::displayObservablesDPDE(vector<Particle> const& /*configuration*/, int iOfStep)
  {
    double totalEnergy = kineticEnergy() + potentialEnergy() + internalEnergy();
    outObservables() << iOfStep * timeStep() 
		    << " " << kineticEnergy()
		    << " " << potentialEnergy()
		    << " " << internalEnergy()
		    << " " << totalEnergy 
		    << std::endl;
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
	outCorrelation() << " " << bendistProfile_(i, midNb) - pow(bendistProfile_.mean(midNb), 2)
			 << std::endl;
      }
  }
 

}
