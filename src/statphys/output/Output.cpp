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
    constTemperature_(input.temperature()),
    constTemperatureLeft_(input.temperature() + input.deltaTemperature()),
    constTemperatureRight_(input.temperature() - input.deltaTemperature()),    
    constDeltaTemperature_(input.deltaTemperature()),
    constGamma_(input.gamma()),
    constXi_(input.xi()),
    constTauBending_(input.tauBending()),
    constExternalForce_(input.externalForce()),
    constNbOfFourier_(input.nbOfFourier()),
    constNbOfHermite_(input.nbOfHermite()),
    totalEnergy_(0),
    totalVirial_(0),
    temperature_(0),
    decorrelationNbOfSteps_(input.decorrelationNbOfSteps()),
    nbOfAutocoPts_(input.nbOfAutocoPts()),
    doFinalFlow_(input.doFinalFlow()),
    doFinalVelocity_(input.doFinalVelocity()),
    obsKineticEnergy_(nullptr),
    obsPotentialEnergy_(nullptr),
    obsPressure_(nullptr),
    obsInternalEnergy_(nullptr),
    obsInternalTemperature_(nullptr),
    obsVelocity_(nullptr),
    obsForce_(nullptr),
    obsLength_(nullptr),
    obsMidFlow_(nullptr),
    obsSumFlow_(nullptr),
    obsModiFlow_(nullptr),
    observables_(),
    kinTempProfile_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfParticles_),
    potTempTopProfile_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfParticles_),
    potTempBotProfile_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfParticles_),
    bendistProfile_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfParticles_),
    flowProfile_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfParticles_),
    cvBasis_(input)
  {    
    if (input.doKineticEnergy()) obsKineticEnergy_ = addObservable(input, "kineticEnergy.txt");
    if (input.doPotentialEnergy()) obsPotentialEnergy_ = addObservable(input, "potentialEnergy.txt");
    if (input.doPressure()) obsPressure_ = addObservable(input, "pressure.txt");
    if (input.doInternalEnergy()) obsInternalEnergy_ = addObservable(input, "internalEnergy.txt");
    if (input.doInternalTemperature()) obsInternalTemperature_ = addObservable(input, "internalTemperature.txt");
    if (input.doVelocity()) obsVelocity_ = addObservable(input, "velocity.txt");
    if (input.doForce()) obsForce_ = addObservable(input, "force.txt");
    if (input.doLength()) obsLength_ = addObservable(input, "length.txt");
    if (input.doMidFlow())obsMidFlow_  = addObservable(input, "midFlow.txt");
    if (input.doSumFlow()) obsSumFlow_ = addObservable(input, "sumFlow.txt");
    if (input.doModiFlow()) obsModiFlow_ = addObservable(input, "modiFlow.txt");
    
    if (input.doOutThermo())
    {
      outThermo_       = std::make_shared<ofstream>(input.outputFolderName() + "thermo.txt");
      outThermo() << "# 1:time  2:position  3:momentum  4:internalEnergy  5:kineticEnergy  6:potentialEnergy  7:totalEnergy  8:kineticTemperature  9:internalTemperature  10:pressure  11:averageRejection  12:negativeEnergiesCount" << endl;
      //outThermo() << "# time kineticEnergy potentialEnergy energy temperature pressure" << endl;
      outMeanThermo_       = std::make_shared<ofstream>(input.outputFolderName() + "meanThermo.txt");
      outMeanThermo() << "# 1:time  2:kineticEnergy  3:potentialEnergy  4:internalEnergy  5:kineticTemperature  6:internalTemperature  7:pressure" << endl;
    }
    if (input.doOutParticles())
    {
      outParticles_ = std::make_shared<ofstream>(input.outputFolderName() + "particles.txt");
      outParticles() << "# time index position momentum kineticEnergy potentialEnergy energy force" << endl;
    }
    if (input.doOutXMakeMol()) outXMakeMol_ = std::make_shared<ofstream>(input.outputFolderName() + "xmakemol.xyz");
    if (input.doOutBackUp()) outBackUp_ = std::make_shared<ofstream>(input.outputFolderName() + "backUp.txt");
    
    if (input.doFinalVelocity()) outFinalVelocity_     = std::make_shared<ofstream>(input.simuTypeName() + "finalVelocity.txt", std::ofstream::app);
    if (input.doFinalFlow()) outFinalFlow_         = std::make_shared<ofstream>(input.simuTypeName() + "finalFlow.txt", std::ofstream::app);
    if (input.doOutChain()) 
    {
      outFinalProfile_      = std::make_shared<ofstream>(input.outputFolderName() + "finalProfile.txt");
      outBeam_              = std::make_shared<ofstream>(input.outputFolderName() + "beamChain.txt");
      outChainVelocities_   = std::make_shared<ofstream>(input.outputFolderName() + "velocitiesChain.txt");
      outChainVelocities() << "# time i=0 i=N/4 i=N/2 i=3N/4 i=N-1" << endl;      
      outProfile_           = std::make_shared<ofstream>(input.outputFolderName() + "profile.txt");
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
    cout << " Observables :";
    for (int iOfObs = 0; iOfObs < nbOfObservables(); iOfObs++)
      cout << "  " << observables(iOfObs)->outPath();
    cout << endl << endl;

  }
  
  Output::~Output()
  {
    for (int iOfObservable = 0; iOfObservable < nbOfObservables(); iOfObservable++)
      delete observables(iOfObservable);
  }
  
  Observable* Output::addObservable(const Input& input, const string& outPath)
  {
    Observable* obsPtr = new Observable(input, outPath);
    observables_.push_back(obsPtr);
    return obsPtr;
  }

  //----- various assessors --------

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
  
  int const& Output::dimension() const
  {return dimension_;}
  
  double const& Output::latticeParameter() const
  {return latticeParameter_;}
  
  int Output::nbOfObservables() const
  {return observables_.size();}
  
  vector<Observable*>& Output::observables()
  {return observables_;}
  
  Observable* Output::observables(int iOfObservable)
  {return observables_[iOfObservable];}
  
  

  const double& Output::kineticEnergy() const
  {return obsKineticEnergy().currentValue();}

  double& Output::kineticEnergy()
  {return obsKineticEnergy().currentValue();}

  const double& Output::potentialEnergy() const
  {return obsPotentialEnergy().currentValue();}

  double& Output::potentialEnergy()
  {return obsPotentialEnergy().currentValue();}

  const double& Output::pressure() const
  {return obsPressure().currentValue();}
  
  double& Output::pressure()
  {return obsPressure().currentValue();}

  const double& Output::totalVirial() const
  {return totalVirial_;}
  
  double& Output::totalVirial()
  {return totalVirial_;}
  
  const double& Output::totalEnergy() const
  {return totalEnergy_;}
  
  double& Output::totalEnergy()
  {return totalEnergy_;}
  
  const double& Output::temperature() const
  {return temperature_;}
  
  double& Output::temperature()
  {return temperature_;}
  
  
  Observable& Output::obsKineticEnergy()
  {return *obsKineticEnergy_;}
  Observable const& Output::obsKineticEnergy() const
  {return *obsKineticEnergy_;}
  Observable& Output::obsPotentialEnergy()
  {return *obsPotentialEnergy_;}
  Observable const& Output::obsPotentialEnergy() const
  {return *obsPotentialEnergy_;}
  Observable& Output::obsPressure()
  {return *obsPressure_;}
  Observable const& Output::obsPressure() const
  {return *obsPressure_;}
  Observable& Output::obsInternalEnergy()
  {return *obsInternalEnergy_;}
  Observable const& Output::obsInternalEnergy() const
  {return *obsInternalEnergy_;}
  Observable& Output::obsInternalTemperature()
  {return *obsInternalTemperature_;}
  Observable const& Output::obsInternalTemperature() const
  {return *obsInternalTemperature_;}
  Observable& Output::obsVelocity()
  {return *obsVelocity_;}
  Observable const& Output::obsVelocity() const
  {return *obsVelocity_;}
  Observable& Output::obsForce()
  {return *obsForce_;}
  Observable const& Output::obsForce() const
  {return *obsForce_;}
  Observable& Output::obsLength()
  {return *obsLength_;}
  Observable const& Output::obsLength() const
  {return *obsLength_;}
  Observable& Output::obsMidFlow()
  {return *obsMidFlow_;}
  Observable const& Output::obsMidFlow() const
  {return *obsMidFlow_;}
  Observable& Output::obsSumFlow()
  {return *obsSumFlow_;}
  Observable const& Output::obsSumFlow() const
  {return *obsSumFlow_;}
  Observable& Output::obsModiFlow()
  {return *obsModiFlow_;}
  Observable const& Output::obsModiFlow() const
  {return *obsModiFlow_;}
    

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
  {return obsInternalEnergy().currentValue();}

  double& Output::internalEnergy()
  {return obsInternalEnergy().currentValue();}

  const double& Output::internalTemperature() const
  {return obsInternalTemperature().currentValue();}

  double& Output::internalTemperature()
  {return obsInternalTemperature().currentValue();}

  const double& Output::rejectionCount() const 
  {return rejectionCount_;}

  double& Output::rejectionCount() 
  {return rejectionCount_;}

  const double& Output::negativeEnergiesCount() const 
  {return negativeEnergiesCount_;}

  double& Output::negativeEnergiesCount() 
  {return negativeEnergiesCount_;}
  
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
    outThermo() << iOfStep * timeStep();
    if (obsKineticEnergy_) outThermo() << " " << kineticEnergy();
    outThermo() << " " << potentialEnergy()
                << " " << totalEnergy()
                << " " << temperature()
                << " " << pressure()
                << std::endl;
  }

  void Output::displayParticles(System const& syst, long int iOfStep)
  {                     
    //for (ParticleIterator it = syst.begin(); !syst.finished(it); ++it)
    for (int iOfParticle = 0; iOfParticle < nbOfParticles_; iOfParticle++)
      outParticles() << iOfStep * timeStep()
                     << " " << iOfParticle
                     << " " << syst(iOfParticle).position()
                     << " " << syst(iOfParticle).momentum()
                     << " " << syst(iOfParticle).kineticEnergy()
                     << " " << syst(iOfParticle).potentialEnergy()
                     << " " << syst(iOfParticle).energy()
                     << " " << syst(iOfParticle).force()
                     << endl; 
  }
  
  ///
  ///-- keep the full current configuration in order to restart from it ---
  void Output::displayBackUp(System const& syst, long int iOfStep)
  {
    outBackUp() << "Time = " << iOfStep * timeStep() << endl;
    int Dim = dimension_;
    for (int iOfParticle = 0; iOfParticle < nbOfParticles_; iOfParticle++)
    {
      for (int dim = 0; dim < Dim; dim++)
        outBackUp() << syst(iOfParticle).position(dim) << " ";
      for (int dim = 0; dim < Dim; dim++)
        outBackUp() << syst(iOfParticle).momentum(dim) << " ";
      outBackUp() << syst(iOfParticle).internalEnergy() << " ";
      outBackUp() << endl;
    }
    outBackUp() << endl;
  }



  void Output::displayThermoVariablesDPDE(System const& syst, long int iOfStep)
  {
    //double totalEnergy = kineticEnergy() + potentialEnergy() + internalEnergy();
    outThermo() << iOfStep * timeStep()
		     << " " << syst(0).position(0) 
		     << " " << syst(0).momentum(0) 
		     << " " << internalEnergy() 
         << " " << kineticEnergy()
         << " " << potentialEnergy()
		     << " " << totalEnergy()
		     << " " << temperature()
		     << " " << 1./internalTemperature()
		     << " " << pressure()
		     << " " << rejectionCount()
		     << " " << negativeEnergiesCount()
		     << std::endl;
    outMeanThermo() << iOfStep * timeStep()
			   << " " << obsKineticEnergy().mean() 
         << " " << obsPotentialEnergy().mean() 
			   << " " << obsInternalEnergy().mean() 
			   << " " << 2*obsKineticEnergy().mean()/(dimension_ * nbOfParticles_)
			   << " " << 1./obsInternalTemperature().mean() 
			   //<< " " << (2*obsKineticEnergy().mean()+obsTotalVirial().mean())/(dimension_*nbOfParticles_*pow(latticeParameter_, dimension_))
         << " " << obsPressure().mean()
			   << std::endl;
  }
  
  // ------------ final outputs ----------------
  
    void Output::displayFinalVelocity()
  {
    outFinalVelocity() << std::left << setw(10) << finalTime()
                        << " " << setw(5) << timeStep()
                        << " " << setw(6) << nbOfParticles()
                        << " " << setw(4) << constTemperature_
                        << " " << setw(6) << constExternalForce_
                        << " " << setw(3) << constNbOfFourier_
                        << " " << setw(3) << constNbOfHermite_
                        << " " << setw(12) << obsVelocity().mean()
                        << " " << setw(12) << obsVelocity().variance()
                        << std::endl;
  }

  void Output::finalDisplayCorrelations()
  {
    //cout << "Velocity : The correlation in 0 is " << floor((2 * obsVelocity().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsVelocity().variance() * 10000)/100 << "% of the variance" << endl;
    //cout << "SumFlow : The correlation in 0 is " << floor((2 * obsSumFlow().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsSumFlow().variance() * 10000)/100 << "% of the variance" << endl;
    //cout << "ModiFlow : The correlation in 0 is " << floor((2 * obsModiFlow().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsModiFlow().variance() * 10000)/100 << "% of the variance" << endl;
        
    for (int iOfObservable=0; iOfObservable < nbOfObservables(); iOfObservable++)
      observables(iOfObservable)->displayCorrelations(nbOfSteps());
  }
  
  bool Output::hasControlVariate() const
  {return cvBasis_.basis_;}


}
