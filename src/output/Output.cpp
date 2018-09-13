#include "simol/output/Output.hpp"

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
    parameters_(input),
    totalVirial_(0),
    temperature_(0),
    decorrelationNbOfSteps_(input.decorrelationNbOfSteps()),
    shortDecorrelationNbOfSteps_(input.shortDecorrelationNbOfSteps()),
    nbOfAutocoPts_(input.nbOfAutocoPts()),
    nbOfShortAutocoPts_(input.nbOfShortAutocoPts()),
    doOutChain_(input.doOutChain()),
    doFinalFlux_(input.doFinalFlux()),
    doFinalLength_(input.doFinalLength()),
    doFinalVelocity_(input.doFinalVelocity()),
    doFinalLagrangeMultiplier_(input.doFinalLagrangeMultiplier()),
    doDPDE_(input.doDPDE()),
    doXMakeMol_(input.doXMakeMol()),
    obsKineticEnergy_(nullptr),
    obsPotentialEnergy_(nullptr),
    obsTotalEnergy_(nullptr),
    obsPressure_(nullptr),
    obsInternalEnergy_(nullptr),
    obsInternalTemperature_(nullptr),
    obsVelocity_(nullptr),
    obsForce_(nullptr),
    obsLength_(nullptr),
    obsMidFlux_(nullptr),
    obsSumFlux_(nullptr),
    obsModiFlux_(nullptr),
    obsLagrangeMultiplier_(nullptr),
    observables_(),
    kinTempProfile_(shortDecorrelationNbOfSteps(), timeStep(), nbOfShortAutocoPts(), nbOfParticles_),
    potTempTopProfile_(shortDecorrelationNbOfSteps(), timeStep(), nbOfShortAutocoPts(), nbOfParticles_),
    potTempBotProfile_(shortDecorrelationNbOfSteps(), timeStep(), nbOfShortAutocoPts(), nbOfParticles_),
    bendistProfile_(shortDecorrelationNbOfSteps(), timeStep(), nbOfShortAutocoPts(), nbOfParticles_),
    fluxProfile_(shortDecorrelationNbOfSteps(), timeStep(), nbOfShortAutocoPts(), nbOfParticles_),
    modiFluxProfile_(shortDecorrelationNbOfSteps(), timeStep(), nbOfShortAutocoPts(), nbOfParticles_),
    extFluxProfile_(shortDecorrelationNbOfSteps(), timeStep(), nbOfShortAutocoPts(), nbOfParticles_)
  {    
    // creates the observables whch are needed, depending on the type of dynamics, system, ...
    // the observable are created with a decorrelation time which is either standard or short 
    // the correlation time is used in the estimation of the asymptotic variance of the observable
    for (int idObs = 0; idObs < nbOfIdObs; idObs++)
      if (idObs == idModiFlux || idObs == idMidFlux || idObs == idLagrangeMultiplier)
        addObservable(input, idObs, shortDecorrelationNbOfSteps(), nbOfShortAutocoPts());
      else
        addObservable(input, idObs, decorrelationNbOfSteps(), nbOfAutocoPts());
    
    // prepares the output file for thermodynamic quantities
    if (input.doOutThermo())
      {
        //-- instantaneous values --
        outThermo_       = std::make_shared<ofstream>(input.outputFolderName() + "thermo.txt");
        if (obsKineticEnergy_)
          {
            if (doDPDE_)
              outThermo() << "# 1:time  2:kineticEnergy  3:potentialEnergy  4:totalEnergy  5:temperature 6:pressure 7:internalEnergy 8:fieldForInternalTemperature 9:negativeEnergiesCountFD 10:negativeEnergiesCountThermal" << endl;
            else
              outThermo() << "# 1:time  2:kineticEnergy  3:potentialEnergy  4:totalEnergy  5:temperature 6:pressure" << endl;
          }
        else
          outThermo() << "# 1:time  2:potentialEnergy  3:totalEnergy  4:temperature 5:pressure" << endl;
        //-- current estimates of averages --
        outMeanThermo_       = std::make_shared<ofstream>(input.outputFolderName() + "meanThermo.txt");
        if (doDPDE_)
          outMeanThermo() << "# 1:time  2:potentialEnergy  3:kineticEnergy  4:temperature 5:pressure 6:internalEnergy 7:internalTemperature 8:averageRejectionRateFD 9:negativeEnergiesCountFD 10:averageRejectionRateThermal 11:negativeEnergiesCountThermal" << endl;
        else
          outMeanThermo() << "# 1:time  2:potentialEnergy  3:kineticEnergy  4:temperature 5:pressure" << endl;
      }

    if (input.doOutParticles())
    {
      outParticles_ = std::make_shared<ofstream>(input.outputFolderName() + "particles.txt");
      outParticles() << "# time index position momentum kineticEnergy potentialEnergy energy force type" << endl;
    }
    if (input.doXMakeMol()) outXMakeMol_        = std::make_shared<ofstream>(input.outputFolderName() + "xmakemol.xyz");
    if (input.doOutBackUp()) outBackUp_            = std::make_shared<ofstream>(input.outputFolderName() + "backUp.txt");
    
    if (doFinalLength()) outFinalLength_ = std::make_shared<ofstream>(input.simuTypeName() + "finalLength.txt", std::ofstream::app);
    if (doFinalVelocity()) outFinalVelocity_ = std::make_shared<ofstream>(input.simuTypeName() + "finalVelocity.txt", std::ofstream::app);
    if (doFinalFlux()) outFinalFlux_         = std::make_shared<ofstream>(input.simuTypeName() + "finalFlux.txt", std::ofstream::app);
    if (doFinalLagrangeMultiplier()) outFinalLagrangeMultiplier_= std::make_shared<ofstream>(input.simuTypeName() + "finalLagrangeMultiplier.txt", std::ofstream::app);
    if (doOutChain()) 
    {
      outFinalProfile_      = std::make_shared<ofstream>(input.outputFolderName() + "finalProfile.txt");
      outInstantProfile_           = std::make_shared<ofstream>(input.outputFolderName() + "instantProfile.txt");
      outInstantProfile() << "# time i=1 i=N/4 i=N/2 i=3N/4 i=N-2" << endl;  
      outProfile_           = std::make_shared<ofstream>(input.outputFolderName() + "profile.txt");
    }
    
    //--------- announce where input/outputs are ----------------------------
    std::cout << endl;
    std::cout << "Output written in " << input.outputFolderName() << std::endl;
    std::cout << "Final output written in " << input.simuTypeName() << std::endl;

    if (outParticles_ && !outParticles_->is_open())
      throw std::runtime_error("The output file does not exist! Please add it manually. ("+input.outputFolderName()+")");

    //-- copy read input into file --
    std::ofstream outInput(input.outputFolderName() + "inputFile.yaml");
    outInput << input.inputFlux().rdbuf();

    //------------------ screen output for control --------------------------
    cout << endl;
    cout << "-----------------------------------------------------" << endl;
    cout << endl;
    cout << " Simulation of : " << endl;
    cout << "    System    : " << input.systemName() << endl;
    cout << "    Dynamics  : " << input.dynamicsName() << endl;
    cout << "    ExternalPotential : " << input.externalPotentialName() << endl;
    cout << "    PairPotential : " << input.pairPotentialName() << endl;
    cout << endl;
    cout << " Number of particles  : " << nbOfParticles_ << endl;
    if (nbOfSteps_ < 1e6)
      cout << " Number of steps      : " << nbOfSteps_ << endl;
    else
      cout << " Number of steps      : " << nbOfSteps_ / 1e6 << "e6" << endl;
    cout << " Time step            : " << timeStep_ << endl;
    cout << " Simulation time      : " << nbOfSteps_*timeStep_ << endl;
    //-- specific screen outputs --
    if (doDPDE())
    {
      cout << endl;
      if (input.MTSfrequency() != 1)
        cout << " -- Multiple timestepping (stoch. part DPDE): " << input.MTSfrequency() << endl;
      if (input.doMetropolis())
        cout << " With Metropolis correction" << endl;
      else 
        cout << " No Metropolis correction" << endl;
      if (input.doProjectionDPDE())
        cout << " With projection for energies" << endl;
      else 
        cout << " No projection" << endl;
      cout << endl;
      if (input.heatCapacityEinstein() > input.heatCapacity())
      {
        cout << " -- blended Einstein model for micro EOS " << endl; 
        cout << " Baseline heat capacity        : " << input.heatCapacity() << endl;
        cout << " Limiting heat capacity        : " << input.heatCapacityEinstein() << endl;
        cout << " Einstein temperature          : " << input.einsteinTemperature() << endl;
      }
      else
      {
        cout << " -- classical micro EOS " << endl;
        cout << " Heat capacity        : " << input.heatCapacity() << endl;
      }
      if ( (input.initialInternalTemperature() != input.temperature()) & (input.thermalizationNbOfSteps() > 0) )
      {
        cout << endl;
        cout << " NONEQUILIBRIUM equilibration dynamics" << endl;
      }
    }
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
  
  ///
  /// memory addresses corresponding to the created observables are made free at the end of the simulation
  Output::~Output()
  {
    for (int iOfObservable = 0; iOfObservable < nbOfObservables(); iOfObservable++)
      delete observables(iOfObservable);
  }
  
  /// Returns the pointer to the observable corresponding to the identifier "idObs"
  /// See Tools.hpp for the definition of these identifiers
  Observable*& Output::getObservablePtr(int idObs)
  {
    switch (idObs)
    {
      case idKineticEnergy : return obsKineticEnergy_;
      case idPotentialEnergy : return obsPotentialEnergy_;
      case idTotalEnergy : return obsTotalEnergy_;
      case idPressure : return obsPressure_;
      case idInternalEnergy : return obsInternalEnergy_;
      case idInternalTemperature : return obsInternalTemperature_;
      case idVelocity : return obsVelocity_;
      case idForce : return obsForce_;
      case idLength : return obsLength_;
      case idMidFlux : return obsMidFlux_;
      case idSumFlux : return obsSumFlux_;
      case idModiFlux : return obsModiFlux_;
      case idLagrangeMultiplier: return obsLagrangeMultiplier_;
      default : throw std::runtime_error("This observable id corresponds to no observable !");
    }
  }
  
  /// Creation of an observable
  /// idObs gives the observable under consideration (see Tools.hpp and getObservablePtr)
  void Output::addObservable(const Input& input, int idObs, int decorrelationNbOfSteps0, int nbOfAutocoPts0)
  {
    Observable* obsPtr = nullptr;
    if (input.doObservable(idObs))
    {
      obsPtr = new Observable(input, idObs, decorrelationNbOfSteps0, nbOfAutocoPts0);
      getObservablePtr(idObs) = obsPtr;
      observables_.push_back(obsPtr);
    }
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

  const double& Output::totalEnergy() const
  {return obsTotalEnergy().currentValue();}

  double& Output::totalEnergy()
  {return obsTotalEnergy().currentValue();}
  
  const double& Output::pressure() const
  {return obsPressure().currentValue();}
  
  double& Output::pressure()
  {return obsPressure().currentValue();}

  const double& Output::totalVirial() const
  {return totalVirial_;}
  
  double& Output::totalVirial()
  {return totalVirial_;}
  
  const double& Output::length() const
  {return obsLength().currentValue();}
  
  double& Output::length()
  {return obsLength().currentValue();}
  
  const double& Output::velocity() const
  {return obsVelocity().currentValue();}
  
  double& Output::velocity()
  {return obsVelocity().currentValue();}
  
  const double& Output::force() const
  {return obsForce().currentValue();}
  
  double& Output::force()
  {return obsForce().currentValue();}
  
  const double& Output::lagrangeMultiplier() const
  {return obsLagrangeMultiplier().currentValue();}
  
  double& Output::lagrangeMultiplier()
  {return obsLagrangeMultiplier().currentValue();}
  
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
  Observable& Output::obsTotalEnergy()
  {return *obsTotalEnergy_;}
  Observable const& Output::obsTotalEnergy() const
  {return *obsTotalEnergy_;}
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
  Observable& Output::obsMidFlux()
  {return *obsMidFlux_;}
  Observable const& Output::obsMidFlux() const
  {return *obsMidFlux_;}
  Observable& Output::obsSumFlux()
  {return *obsSumFlux_;}
  Observable const& Output::obsSumFlux() const
  {return *obsSumFlux_;}
  Observable& Output::obsModiFlux()
  {return *obsModiFlux_;}
  Observable const& Output::obsModiFlux() const
  {return *obsModiFlux_;}
  Observable& Output::obsLagrangeMultiplier()
  {return *obsLagrangeMultiplier_;}
  Observable const& Output::obsLagrangeMultiplier() const
  {return *obsLagrangeMultiplier_;}
    

  const double& Output::timeStep() const
  {return timeStep_;}

  double& Output::timeStep()
  {return timeStep_;}

  /*int & Output::decorrelationNbOfSteps()
  {return decorrelationNbOfSteps_;}*/

  const int& Output::decorrelationNbOfSteps() const
  {return decorrelationNbOfSteps_;}

  double Output::decorrelationTime() const
  {return decorrelationNbOfSteps() * timeStep();}
  
  const int& Output::shortDecorrelationNbOfSteps() const
  {return shortDecorrelationNbOfSteps_;}

  double Output::shortDecorrelationTime() const
  {return shortDecorrelationNbOfSteps() * timeStep();}

  const int& Output::nbOfAutocoPts() const
  {return nbOfAutocoPts_;}

  double Output::autocoPtsPeriod() const
  {return decorrelationTime() / nbOfAutocoPts();}
  
  const int& Output::nbOfShortAutocoPts() const
  {return nbOfShortAutocoPts_;}

  double Output::shortAutocoPtsPeriod() const
  {return shortDecorrelationTime() / nbOfAutocoPts();}

  //---- specific for DPDE -----
  
  const double& Output::internalEnergy() const
  {return obsInternalEnergy().currentValue();}

  double& Output::internalEnergy()
  {return obsInternalEnergy().currentValue();}

  const double& Output::internalTemperature() const
  {return obsInternalTemperature().currentValue();}

  double& Output::internalTemperature()
  {return obsInternalTemperature().currentValue();}

  const double& Output::rejectionCountFD() const 
  {return rejectionCountFD_;}

  double& Output::rejectionCountFD() 
  {return rejectionCountFD_;}

  const double& Output::negativeEnergiesCountFD() const 
  {return negativeEnergiesCountFD_;}

  double& Output::negativeEnergiesCountFD() 
  {return negativeEnergiesCountFD_;}

  const double& Output::rejectionCountThermal() const 
  {return rejectionCountThermal_;}

  double& Output::rejectionCountThermal() 
  {return rejectionCountThermal_;}

  const double& Output::negativeEnergiesCountThermal() const 
  {return negativeEnergiesCountThermal_;}

  double& Output::negativeEnergiesCountThermal() 
  {return negativeEnergiesCountThermal_;}
  
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

  ///
  /// writes some relevent thermodynamic quantities in the output file "thermo.txt"
  void Output::displayThermoVariables(long int iOfStep)
  {
    //---- instantaneous values ----
    outThermo() << iOfStep * timeStep();
    if (obsKineticEnergy_) outThermo() << " " << kineticEnergy();
    outThermo() << " " << potentialEnergy()
                << " " << totalEnergy()
                << " " << temperature()
                << " " << pressure();
    if (doDPDE_) // add fields for DPDE
    {
      outThermo() << " " << internalEnergy() 
		  << " " << 1./internalTemperature()  // in fact, field used to compute internal temp. as ratio of averages when general microEOS; otherwise instantaneous inverse internal temperature (should not be averaged as such...)
		  << " " << negativeEnergiesCountFD()
		  << " " << negativeEnergiesCountThermal();
    }
    outThermo() << std::endl;
    
    //---- current estimates of the averages ----
    outMeanThermo() << iOfStep * timeStep();
    if (obsKineticEnergy_) outMeanThermo() << " " << obsKineticEnergy().mean(); 
    outMeanThermo() << " " << obsPotentialEnergy().mean(); 
    if (obsKineticEnergy_) outMeanThermo() << " " << 2*obsKineticEnergy().mean()/(dimension_ * nbOfParticles_);
    outMeanThermo() << " " << obsPressure().mean();
    if (doDPDE_) // add fields for DPDE
    {
      outMeanThermo() << " " << obsInternalEnergy().mean() 
		      << " " << 1./obsInternalTemperature().mean() // for general microEOS, use: obsInternalEnergy().mean()/nbOfParticles_/(1+obsInternalTemperature().mean())
		      << " " << rejectionCountFD()
		      << " " << negativeEnergiesCountFD()
		      << " " << rejectionCountThermal()
		      << " " << negativeEnergiesCountThermal() ;
    }
    outMeanThermo() << std::endl;
    
  }

  ///
  /// writes one line of output per particle into the ouput file "particle.txt"
  void Output::displayParticles(System const& syst, long int iOfStep)
  {
    for (int iOfParticle = 0; iOfParticle < nbOfParticles_; iOfParticle++)
      outParticles() << iOfStep * timeStep()
                     << " " << iOfParticle
                     << " " << syst.periodicImage(syst(iOfParticle).position()).transpose()
                     << " " << syst(iOfParticle).momentum().transpose()
                     << " " << syst(iOfParticle).kineticEnergy()
                     << " " << syst(iOfParticle).potentialEnergy()
                     << " " << syst(iOfParticle).energy()
                     << " " << syst(iOfParticle).force().transpose()
                     << " " << syst(iOfParticle).type()
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
      if (doDPDE_) outBackUp() << syst(iOfParticle).internalEnergy() << " ";
      outBackUp() << endl;
    }
    outBackUp() << endl;
  }
  
  // ------------ final outputs ----------------
  
  ///
  /// Rarely used, writes statistics relative to a relevent length (or position) in the ouput file "length.txt"
  void Output::displayFinalLength()
  {
    outFinalLength() << std::left << setw(10) << finalTime()
                        << " " << setw(5) << timeStep()
                        << " " << setw(6) << nbOfParticles()
                        << " " << setw(4) << parameters_.temperature()
                        << " " << setw(6) << parameters_.eta()
                        << " " << setw(6) << parameters_.interactionRatio()
                        //<< " " << setw(3) << parameters_.nbOfQModes()
                        //<< " " << setw(3) << parameters_.nbOfPModes()
                        << " " << setw(12);
    
    obsLength().displayFinalValues(outFinalLength());
  }
  
  /// Called at the end of the simulation
  /// Writes a single line with the parameters of the dynamics and estimated statstics of the velocity (or a relevent velocity)
  void Output::displayFinalVelocity()
  {
    outFinalVelocity() << std::left << setw(10) << finalTime()
                        << " " << setw(5) << timeStep()
                        << " " << setw(4) << nbOfParticles()
                        << " " << setw(4) << parameters_.temperature()
                        << " " << setw(8) << parameters_.gamma()
                        << " " << setw(8) << parameters_.eta()
                        << " " << setw(8) << parameters_.coupling()
                        //<< " " << setw(3) << parameters_.nbOfQModes()
                        //<< " " << setw(3) << parameters_.nbOfPModes()
                        << " " << setw(12);
                        
    obsVelocity().displayFinalValues(outFinalVelocity());
  }
  
  /// Called at the end of the simulation
  /// Writes a single line with the parameters of the dynamics and estimated statstics of the Lagrange multiplier
  /// This output is such that the mobility can be estimated as the unbiased drift divided by the mean Lagrange multiplier
  void Output::displayFinalLagrangeMultiplier()
  {    
    cout << "outFinalLagrangeMultiplier_" << std::left << " " << setw(10) << finalTime()
                    << " " << setw(5) << timeStep()
                    << " " << setw(6) << nbOfParticles()
                    << " " << setw(4) << parameters_.temperature()
                    << " " << setw(6) << parameters_.drift()
                    //<< " " << setw(9) << unbiasedDrift
                    << " " << setw(12) << obsLagrangeMultiplier().mean()
                    << " " << setw(12) << obsLagrangeMultiplier().asymptoticVariance()
                    << " " << setw(12) << obsLagrangeMultiplier().asyvarOfAsyvar();
    
    outFinalLagrangeMultiplier() << std::left << setw(10) << finalTime()
                        << " " << setw(5) << timeStep()
                        << " " << setw(6) << nbOfParticles()
                        << " " << setw(4) << parameters_.temperature()
                        << " " << setw(6) << parameters_.drift()
                        //<< " " << setw(9) << unbiasedDrift
                        << " " << setw(12);    
    obsLagrangeMultiplier().displayFinalValues(outFinalLagrangeMultiplier());  
  }

  ///
  /// Writes all the correlation profiles which have been estimated in their respective output files "correlation....txt"
  void Output::finalDisplayCorrelations()
  {   
    for (int iOfObservable=0; iOfObservable < nbOfObservables(); iOfObservable++)
      observables(iOfObservable)->displayCorrelations(nbOfSteps());
  }
  


}
