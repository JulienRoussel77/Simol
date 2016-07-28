#include "simol/statphys/output/Output.hpp"

using std::cout;
using std::endl;
using std::vector;

namespace simol
{

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
    //energyMidFlow_(0),
    //energySumFlow_(0),
    internalTemperature_(0),
    decorrelationNbOfSteps_(input.decorrelationNbOfSteps()),
    nbOfAutocoPts_(input.nbOfAutocoPts()),
    doFinalFlow_(input.doFinalFlow()),
    doFinalVelocity_(input.doFinalVelocity()),
    obsVelocity_(nullptr),
    obsForce_(nullptr),
    obsLength_(nullptr),
    obsMidFlow_(nullptr),
    obsSumFlow_(nullptr),
    observables_(),
    /*velocityCV_(nullptr),
    obsForce_(nullptr),
    obsLength_(nullptr),
    obsMidFlow_(nullptr),
    obsSumFlow_(nullptr),*/
    kinTempProfile_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfParticles_),
    potTempTopProfile_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfParticles_),
    potTempBotProfile_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfParticles_),
    bendistProfile_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfParticles_),
    flowProfile_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfParticles_),
    averageKineticEnergy_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts()),
    averagePotentialEnergy_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts()),
    averageInternalEnergy_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts()),
    averageInternalTemperature_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts()),
    averagePressure_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts())
    
  {
    
    //-- std observables in this file --
    outThermo_       = std::make_shared<ofstream>(input.outputFolderName() + "thermo.txt");
    outThermo() << "# time kineticEnergy potentialEnergy energy temperature pressure" << endl;

    //-- override the standrd observable files for DPDE --
    if (input.dynamicsName() == "DPDE") 
    {
      outThermo_       = std::make_shared<ofstream>(input.outputFolderName() + "thermo.txt");
      outThermo() << "# 1:time  2:position  3:momentum  4:internalEnergy  5:kineticEnergy  6:potentialEnergy  7:totalEnergy  8:kineticTemperature  9:internalTemperature  10:pressure  11:averageRejection  12:negativeEnergiesCount" << endl;
      meanValueObservables_       = std::make_shared<ofstream>(input.outputFolderName() + "mean_thermo.txt");

      meanValueObservables() << "# 1:time  2:kineticEnergy  3:potentialEnergy  4:internalEnergy  5:kineticTemperature  6:internalTemperature  7:pressure" << endl;
    }

    //-- longer outputs if required, e.g. configuration of the system --
    if (printLongPeriodNbOfSteps_ > 0)
    {
      if (input.systemName() == "NBody")
      {
        outParticlesXMakeMol_ = std::make_shared<ofstream>(input.outputFolderName() + "xmakemol.xyz");
        outParticlesFullConfiguration_ = std::make_shared<ofstream>(input.outputFolderName() + "FullConfiguration.txt");
      }
      else
      {
        outParticles_ = std::make_shared<ofstream>(input.outputFolderName() + "particles.txt");
        outParticles() << "# time index position momentum kineticEnergy potentialEnergy energy force" << endl;
      }
    }

    //-- for autocorrelations --
    if (decorrelationNbOfSteps_ > 0)
      outCorrelation_ = std::make_shared<ofstream>(input.outputFolderName() + "correlations.txt");

    //-- average velocities for 1D nonequilibrium --
    if ( (input.systemName() == "Isolated") && (input.dynamicsName() == "Langevin") )
    {
      //outVelocitiesCV_      = std::make_shared<ofstream>(input.outputFolderName() + "velocities.txt");
      //outVelocitiesCV() << "# time b <b> <b2> D <D> <D2> observable <observable> <observable2> LPhi <LPhi> <LPhi2>" << endl;

        obsVelocity_ = addObservable(input, "velocity.txt");
    }

    //-- average forces for 1D nonequilibrium overdamped --
    if ( (input.systemName() == "Isolated") && (input.dynamicsName() == "Overdamped") )
    {
      //outForcesCV_ = std::make_shared<ofstream>(input.outputFolderName() + "forces.txt");
      //outForcesCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
      obsForce_ = addObservable(input, "force.txt");
    }

    //-- used for control variates --
    if (input.controlVariateName() != "None")
      outVelocitiesGenerator_ = std::make_shared<ofstream>(input.outputFolderName() + "velocitiesGenerator.txt");

    //--- to have a summary of the fina, estimated average velocity ---
    if (doFinalVelocity_)
      outFinalVelocity_     = std::make_shared<ofstream>(input.simuTypeName() + "finalVelocity.txt", std::ofstream::app);

    //-- outputs specific to chains --
    if (input.systemName() == "BiChain" || input.systemName() == "TriChain")
    {
      outFinalProfile_      = std::make_shared<ofstream>(input.outputFolderName() + "finalProfile.txt");

      outBeam_              = std::make_shared<ofstream>(input.outputFolderName() + "beamChain.txt");
      
      //outLengthsCV_         = std::make_shared<ofstream>(input.outputFolderName() + "lengths.txt");
      //outLengthsCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
      obsLength_ = addObservable(input, "length.txt");

      outChainVelocities_   = std::make_shared<ofstream>(input.outputFolderName() + "velocitiesChain.txt");
      outChainVelocities() << "# time i=0 i=N/4 i=N/2 i=3N/4 i=N-1" << endl;

      /*outMidFlowCV_         = std::make_shared<ofstream>(input.outputFolderName() + "midFlow.txt");
      outMidFlowCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
      outMidFlowPT_         = std::make_shared<ofstream>(input.outputFolderName() + "midFlowPost.txt");*/
      
      obsMidFlow_ = addObservable(input, "midFlow.txt");

      /*outSumFlowCV_         = std::make_shared<ofstream>(input.outputFolderName() + "sumFlow.txt");
      outSumFlowCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;
      outSumFlowPT_         = std::make_shared<ofstream>(input.outputFolderName() + "sumFlowPost.txt");*/
      
      obsSumFlow_ = addObservable(input, "sumFlow.txt");

      /*outModiFlowCV_         = std::make_shared<ofstream>(input.outputFolderName() + "modiFlow.txt");
      outModiFlowCV() << "# time b <b> <b2> D <D> >D2> observable <observable> >observable2> LPhi <LPhi> <LPhi2>" << endl;*/
      
      obsModiFlow_ = addObservable(input, "modiFlow.txt");
      
      outProfile_           = std::make_shared<ofstream>(input.outputFolderName() + "profile.txt");
      
      outTest_           = std::make_shared<ofstream>(input.outputFolderName() + "test.txt");

      if (doFinalFlow_)
        outFinalFlow_         = std::make_shared<ofstream>(input.simuTypeName() + "finalFlow.txt", std::ofstream::app);
    }

    //--------- announce where input/outputs are ----------------------------
    std::cout << endl;
    std::cout << "Output written in " << input.outputFolderName() << std::endl;
    std::cout << "Final output written in " << input.simuTypeName() << std::endl;

    if (outParticles_ && !outParticles_->is_open())
      throw std::runtime_error("The output file does not exist ! Please add it manually.");

    //-- copy read input into file --
    std::ofstream outInput(input.outputFolderName() + "inputFile.txt");
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
    if (nbOfSteps_ < 1e6)
      cout << " Number of steps      : " << nbOfSteps_ << endl;
    else
      cout << " Number of steps      : " << nbOfSteps_ / 1e6 << "e6" << endl;
    cout << " Time step            : " << timeStep_ << endl;
    //-- specific screen outputs --
    if (input.dynamicsName() == "DPDE")
      cout << " Heat capacity        : " << input.heatCapacity() << endl;
    if (input.systemName() == "NBody" && input.doCellMethod())
    {
      cout << endl;
      cout << " Using the method of cells " << endl;
    }
    cout << " Period between light outputs    : " << input.printPeriodTime() << endl;
    cout << " Period between heavy outputs    : " << input.printLongPeriodTime() << endl;
    cout << endl;
    cout << "-----------------------------------------------------" << endl << endl;

  }
  
  Observable* Output::addObservable(const Input& input, const string& outPath)
  {
    Observable* obsPtr = new Observable(input, outPath);
    observables_.push_back(obsPtr);
    return obsPtr;
  }

  //----- various assessors --------

  ofstream & Output::outThermo()
  {return *outThermo_;}

  ofstream & Output::outParticles()
  {return *outParticles_;}

  ofstream & Output::outFinalVelocity()
  {return *outFinalVelocity_;}

  ofstream & Output::outCorrelation()
  {return *outCorrelation_;}

  ofstream & Output::meanValueObservables()
  {return *meanValueObservables_;}

  //ofstream & Output::outVelocitiesCV()
  //{return *outVelocitiesCV_;}

  ofstream & Output::outVelocitiesGenerator()
  {return *outVelocitiesGenerator_;}

  /*ofstream & Output::outForcesCV()
  {return *outForcesCV_;}

  ofstream & Output::outLengthsCV()
  {return *outLengthsCV_;}*/

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

  const long int& Output::nbOfSteps() const
  {return nbOfSteps_;}

  double Output::finalTime() const
  {return nbOfSteps_ * timeStep_;}
  
  int Output::nbOfObservables() const
  {return observables_.size();}
  
  Observable* Output::observables(int iOfObservable)
  {return observables_[iOfObservable];}
  
  

  const double& Output::kineticEnergy() const
  {return kineticEnergy_;}

  double& Output::kineticEnergy()
  {return kineticEnergy_;}

  const double& Output::potentialEnergy() const
  {return potentialEnergy_;}

  double& Output::potentialEnergy()
  {return potentialEnergy_;}

  const double& Output::totalVirial() const
  {return totalVirial_;}

  double& Output::totalVirial()
  {return totalVirial_;}

  double Output::energy() const
  {return kineticEnergy_ + potentialEnergy_;}

  double Output::temperature() const
  {return 2 * kineticEnergy_ / (dimension_ * nbOfParticles_); }

  double Output::pressure() const
  {return (2 * kineticEnergy_ + totalVirial_) / (dimension_ * nbOfParticles_ * pow(latticeParameter_, dimension_));}

  bool Output::doComputeCorrelations() const
  {return decorrelationNbOfSteps() > 0;}

  Observable& Output::obsVelocity()
  {return *obsVelocity_;}

  Observable& Output::obsForce()
  {return *obsForce_;}

  Observable& Output::obsLength()
  {return *obsLength_;}

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

  //---- specific for DPDE -----
  
  const double& Output::internalEnergy() const
  {return internalEnergy_;}

  double& Output::internalEnergy()
  {return internalEnergy_;}

  const double& Output::rejectionCount() const 
  {return rejectionCount_;}

  double& Output::rejectionCount() 
  {return rejectionCount_;}

  const double& Output::negativeEnergiesCount() const 
  {return negativeEnergiesCount_;}

  double& Output::negativeEnergiesCount() 
  {return negativeEnergiesCount_;}

  const double& Output::internalTemperature() const
  {return internalTemperature_;}

  double& Output::internalTemperature()
  {return internalTemperature_;}

  //---------------- check whether outputs should be performed --------------------

  bool Output::doOutput(long int iOfStep) const
  {
    return (printPeriodNbOfSteps() > 0 && iOfStep % printPeriodNbOfSteps() == 0);
  }

  bool Output::doLongPeriodOutput(long int iOfStep) const
  {
    return (printLongPeriodNbOfSteps() > 0 && iOfStep % printLongPeriodNbOfSteps() == 0);
  }

  //----------------- display functions --------------------------

  void Output::displayThermoVariables(long int iOfStep)
  {
    outThermo() << iOfStep * timeStep()
                     << " " << kineticEnergy()
                     << " " << potentialEnergy()
                     << " " << energy()
                     << " " << temperature()
                     << " " << pressure()
                     << std::endl;
  }

  void Output::displayParticles(vector<Particle> const& configuration, long int iOfStep)
  {
    for (int iOfParticle = 0; iOfParticle < nbOfParticles_; iOfParticle++)
      outParticles() << iOfStep * timeStep()
                     << " " << iOfParticle
                     << " " << configuration[iOfParticle].position()
                     << " " << configuration[iOfParticle].momentum()
                     << " " << configuration[iOfParticle].kineticEnergy()
                     << " " << configuration[iOfParticle].potentialEnergy()
                     << " " << configuration[iOfParticle].energy()
                     << " " << configuration[iOfParticle].force()
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
                         << " " << setw(12) << obsVelocity().mean()
                         << " " << setw(12) << obsVelocity().variance()
                         << std::endl;
    }
  }
  
  void Output::appendKineticEnergy(double value, long int iOfStep)
  {
    averageKineticEnergy_.append(value, iOfStep);
  }

  void Output::appendPotentialEnergy(double value, long int iOfStep)
  {
    averagePotentialEnergy_.append(value, iOfStep);
  }
  
  void Output::appendInternalEnergy(double value, long int iOfStep)
  {
    averageInternalEnergy_.append(value, iOfStep);
  }

  void Output::appendInternalTemperature(double value, long int iOfStep)
  {
    averageInternalTemperature_.append(value, iOfStep);
  }

  void Output::appendPressure(double value, long int iOfStep)
  {
    averagePressure_.append(value, iOfStep);
  }

  void Output::displayThermoVariablesDPDE(vector<Particle> const& configuration, long int iOfStep)
  {
    double totalEnergy = kineticEnergy() + potentialEnergy() + internalEnergy();
    outThermo() << iOfStep * timeStep()
		     << " " << configuration[0].position(0) 
		     << " " << configuration[0].momentum(0) 
		     << " " << internalEnergy() 
                     << " " << kineticEnergy()
                     << " " << potentialEnergy()
		     << " " << totalEnergy
		     << " " << temperature()
		     << " " << 1./internalTemperature()
		     << " " << pressure()
		     << " " << rejectionCount()
		     << " " << negativeEnergiesCount()
		     << std::endl;
    meanValueObservables() << iOfStep * timeStep()
			   << " " << averageKineticEnergy_.mean() 
      			   << " " << averagePotentialEnergy_.mean() 
			   << " " << averageInternalEnergy_.mean() 
			   << " " << 2*averageKineticEnergy_.mean()/(dimension_ * nbOfParticles_)
			   << " " << 1./averageInternalTemperature_.mean() 
			   << " " << (2*averageKineticEnergy_.mean()+averagePressure_.mean())/(dimension_*nbOfParticles_*pow(latticeParameter_, dimension_))
			   << std::endl;
  }

  //------------ autocorrelations ------------------
  
  void Output::finalDisplayCorrelationsDPDE()
  {
    for (int iOfSpan = 0; iOfSpan < nbOfAutocoPts(); iOfSpan++)
      {
        //-- autocoPtsPeriod(): time between successive correlation values; may be different from the timestep if some subsampling is specified
        outCorrelation() << iOfSpan * autocoPtsPeriod() 
            << " " << averageKineticEnergy_.unbiasedCorrelationAtSpan(iOfSpan)
            << " " << averagePotentialEnergy_.unbiasedCorrelationAtSpan(iOfSpan)
            << " " << averageInternalEnergy_.unbiasedCorrelationAtSpan(iOfSpan)
            << " " << averageInternalTemperature_.unbiasedCorrelationAtSpan(iOfSpan)
            << endl;
      }
  }

  void Output::finalDisplayCorrelations()
  {
    //cout << "Velocity : The correlation in 0 is " << floor((2 * obsVelocity().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsVelocity().variance() * 10000)/100 << "% of the variance" << endl;
    //cout << "SumFlow : The correlation in 0 is " << floor((2 * obsSumFlow().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsSumFlow().variance() * 10000)/100 << "% of the variance" << endl;
    //cout << "ModiFlow : The correlation in 0 is " << floor((2 * obsModiFlow().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsModiFlow().variance() * 10000)/100 << "% of the variance" << endl;
        
    if (doComputeCorrelations())
      for (int iOfObservable=0; iOfObservable < nbOfObservables(); iOfObservable++)
        observables(iOfObservable)->displayCorrelations(nbOfSteps());
    
    /*int midNb = nbOfParticles_ / 2;
    for (int iOfSpan = 0; iOfSpan < nbOfAutocoPts(); iOfSpan++)
    {
      outCorrelation() << iOfSpan * autocoPtsPeriod()
                       << " " << obsVelocity().unbiasedCorrelationAtSpan(iOfSpan)
                       << " " << sqrt(obsVelocity().varCorrelationAtSpan(iOfSpan) / finalTime())
                       << " " << obsForce_->unbiasedCorrelationAtSpan(iOfSpan)
                       << " " << sqrt(obsForce_->varCorrelationAtSpan(iOfSpan) / finalTime())
                       << " " << obsLength_->unbiasedCorrelationAtSpan(iOfSpan)
                       << " " << sqrt(obsLength_->varCorrelationAtSpan(iOfSpan) / finalTime())
                       << " " << obsMidFlow_->unbiasedCorrelationAtSpan(iOfSpan)
                       << " " << sqrt(obsMidFlow_->varCorrelationAtSpan(iOfSpan) / finalTime())
                       << " " << obsSumFlow_->unbiasedCorrelationAtSpan(iOfSpan)
                       << " " << sqrt(obsSumFlow_->varCorrelationAtSpan(iOfSpan) / finalTime())  //#11
                       << " " << obsModiFlow_->unbiasedCorrelationAtSpan(iOfSpan)
                       << " " << sqrt(obsModiFlow_->varCorrelationAtSpan(iOfSpan) / finalTime())
                       << " " << bendistProfile_.unbiasedCorrelationAtSpan(iOfSpan, midNb)
                       << std::endl;
    }*/
  }


}
