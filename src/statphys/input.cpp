#include "input.hpp"

using std::cout; 
using std::endl; 
using std::string;
using std::max;
using std::min;

namespace simol {
	
  ///Transforms a double to a nice string
  /// 1.30000 -> 1.3     1.000 -> 1  
  string doubleToString(double x) 
	{
		string s = to_string (x);
		if (s.find(".") != std::string::npos)
			s.erase ( s.find_last_not_of('0') + 1, string::npos );
		s.erase ( s.find_last_not_of('.') + 1, string::npos );
		return s;
	}
	
	///
	///Reads the next element in the file read
	double readItem(ifstream& in)
	{
		double item;
		in >> item;
		return item;
	}
  
	///Reads the input file using the YAML library
	///Reads the initial conditions in the "settings" file, if indicated
  Input::Input(CommandLine cmd):
		data(YAML::LoadFile(cmd.inputFileName())),
		inputPath_(cmd.inputFileName()),
		inputFlux_(inputPath()),
		inputSettings_(settingsPath())
  {
		assert(data["Physics"]);
    
    if (data["Physics"]["System"]["Position"])
		{
      if (data["Physics"]["System"]["Position"].size() == 2)
      {
				cout << "input is double !" << endl;
				exit(1);
				positionMin_ = data["Physics"]["System"]["Position"][0].as<double>();
				positionMax_ = data["Physics"]["System"]["Position"][1].as<double>();
      }
		}
      
    
		else if (doFileSettings())
		{
			assert(inputFlux_.is_open());
			cout << "Reading the settings from " << settingsPath() << "...";
			initialPositions_ = vector<dvec>(nbOfParticles());
			initialMomenta_ = vector<dvec>(nbOfParticles());
			for (int iOfParticle=0; iOfParticle < (int)nbOfParticles(); iOfParticle++)
			{
				readItem(inputSettings_);
				assert( (int) readItem(inputSettings_) == (int) iOfParticle);
				for (int i=0; i < dimension(); i++)
					initialPositions_[iOfParticle](i) = readItem(inputSettings_);
				for (int i=0; i < dimension(); i++)
					initialMomenta_[iOfParticle](i) = readItem(inputSettings_);
				for (int i=0; i<4; i++)
					readItem(inputSettings_);
				//cout << iOfParticle << " " << initialPositions_[iOfParticle] << " " << initialMomenta_[iOfParticle] << endl;
			}
			cout <<"OK !" << endl;
		}
  }
  
  const std::string& Input::inputPath() const
  {
		return inputPath_;
	}
	
	const std::ifstream& Input::inputFlux() const
	{
		return inputFlux_;
	}
	
	string Input::simuTypeName() const {
    return "../output/"+dynamicsName()+"/"+systemName()+"/"+potentialName()+"/";
  }
	
	string Input::parametersName() const
  {
		string name = simuTypeName();
		
		if (controlVariateName() != "None")
      name += controlVariateName()+"/";
		
		if (data["Output"]["ParametersName"])
			return name + data["Output"]["ParametersName"].as<string>()+"/";
		
		if (dynamicsName() == "BoundaryLangevin")
			name += "N" + to_string(nbOfParticles()) + "_";
		name += "dt" + doubleToString(timeStep()) + "_eta" + doubleToString(eta());
		if (dynamicsName() == "BoundaryLangevin")
			name += "_xi" + doubleToString(xi());
		name += "/";
		return name;
	}

  string Input::outputFolderName() const {
    
		string name = parametersName();		

    if (data["Output"]["FolderName"])
      name += data["Output"]["FolderName"].as<string>()+"/";
		
    return name;
  }
	
	//### Geometry ###

  int Input::dimension() const {return data["Geometry"]["Dimension"].as<int>();}

  double Input::length() const 
  {
		if (data["Geometry"]["Length"])
			return data["Geometry"]["Length"].as<double>();
		else
			return 2 * M_PI;
	}
	
	// ### Mesh/Time ###
  
  double Input::timeStepMin() const 
  {
    if (data["Mesh"]["Time"]["Step"].size() == 2)
      return data["Mesh"]["Time"]["Step"][0].as<double>();
    else
      return data["Mesh"]["Time"]["Step"].as<double>(); 
  }
  
    double Input::timeStepMax() const 
  {
    if (data["Mesh"]["Time"]["Step"].size() == 2)
      return data["Mesh"]["Time"]["Step"][1].as<double>();
    else
      return data["Mesh"]["Time"]["Step"].as<double>(); 
  }

  double Input::timeStep(size_t iOfReplica) const 
  {
    return timeStepMin() * pow(timeStepMax() / timeStepMin(), iOfReplica / max(1., nbOfReplicas()-1.));
  }
  
  ///
	///Number of iterations in the "output" part of the simulation
  size_t Input::nbOfIterations(size_t iOfReplica) const 
  {
    if (data["Mesh"]["Time"]["Number"])
      return data["Mesh"]["Time"]["Number"].as<size_t>();
    else if (data["Mesh"]["Time"]["FinalTime"])
      return data["Mesh"]["Time"]["FinalTime"].as<double>() / timeStep(iOfReplica);
    else
			{cout << "Nb of Iterations not specified !" << endl;exit(1);}
  }
  
  ///
  ///Number of iterations in the thermalization part (temperatures imposed everywhere)
  size_t Input::nbOfThermalIterations(size_t iOfReplica) const 
  {
    if (data["Mesh"]["Time"]["ThermalNumber"])
      return data["Mesh"]["Time"]["ThermalNumber"].as<size_t>();
    else if (data["Mesh"]["Time"]["ThermalTime"])
      return data["Mesh"]["Time"]["ThermalTime"].as<double>() / timeStep(iOfReplica);
    else
			return 0;
  }
  
  ///
  ///Number of iterations in the burning part (no output)
  size_t Input::nbOfBurningIterations(size_t iOfReplica) const 
  {
    if (data["Mesh"]["Time"]["BurningNumber"])
      return data["Mesh"]["Time"]["BurningNumber"].as<size_t>();
    else if (data["Mesh"]["Time"]["BurningTime"])
      return data["Mesh"]["Time"]["BurningTime"].as<double>() / timeStep(iOfReplica);
    else
			return 0;
  }
  
  // ### Mesh/Replica ###
  
  size_t Input::nbOfReplicas() const {
    if (data["Mesh"]["Replicas"]["Number"])
      return data["Mesh"]["Replicas"]["Number"].as<int>();
    else return 1;
  }
  
  // ### Physics/System ###

  string Input::systemName() const {return data["Physics"]["System"]["Name"].as<string>();}
  
  size_t Input::nbOfParticles() const {
    if (data["Physics"]["System"]["Number"])
      return data["Physics"]["System"]["Number"].as<size_t>();
    else return 1;
  }
  
  double Input::mass() const {
    if (data["Physics"]["System"]["Mass"])
      return data["Physics"]["System"]["Mass"].as<double>();
    else return 1;
  } 	
  
  ///
  ///Returns True if the initial confitions are read from the "settings" file
  bool Input::doFileSettings() const
  {
		if (data["Physics"]["System"]["Settings"])
			if (data["Physics"]["System"]["Settings"].as<string>() == "FileSettings")
				return true;
		return false;
	}
	
  string Input::settingsPath() const
  {
		if (data["Physics"]["System"]["SettingsPath"])
			return simuTypeName()+data["Physics"]["System"]["SettingsPath"].as<string>();
		else return parametersName()+"settings/settings";
	}
  
  dvec Input::initialPosition(int const& iOfParticle) const {
    if (data["Physics"]["System"]["Position"])
		{
      if (data["Physics"]["System"]["Position"].size() == 1)
				return data["Physics"]["System"]["Position"].as<double>();
      else
				return positionMin_ + (iOfParticle + .5)/nbOfParticles() * (positionMax_ - positionMin_);
		}
		else if (doFileSettings())
		{
			cout << "using settings for q : " << iOfParticle << "->" << initialPositions_[iOfParticle] << endl;
			return initialPositions_[iOfParticle];
		}
			
		else return dvec(dimension(), 0);
  }   
  
  dvec Input::initialMomentum(int const& iOfParticle) const 
  {
    if (data["Physics"]["System"]["Momentum"])
      return data["Physics"]["System"]["Momentum"].as<double>();
		else if (doFileSettings())
		{
			cout << "using settings for p : " << iOfParticle << "->" << initialMomenta_[iOfParticle] << endl;
			return initialMomenta_[iOfParticle];
		}
    else return dvec(dimension(), 0);
    //else return (i < nbOfParticles()/2)?.2:-.2;
  }  
  
	// ### Physics/Potential ###
  
  string Input::potentialName() const {return data["Physics"]["Potential"]["Name"].as<string>();}

  //Sinusoidal
  double Input::amplitude() const 
  {
    if (data["Physics"]["Potential"]["Amplitude"])
      return data["Physics"]["Potential"]["Amplitude"].as<double>();
    else
      return 1;
  }
  
  //DoubleWell
  double Input::height() const {return data["Physics"]["Potential"]["Height"].as<double>();}
  double Input::interWell() const {return data["Physics"]["Potential"]["InterWell"].as<double>();}
  
  //Harmonic
  double Input::potentialStiffness() const {
    if (data["Physics"]["Potential"]["Stiffness"])
      return data["Physics"]["Potential"]["Stiffness"].as<double>();
    else
      return 1;
  }
  
  //Quadratic
  double Input::potentialAlpha() const {
    if (data["Physics"]["Potential"]["Alpha"])
      return data["Physics"]["Potential"]["Alpha"].as<double>();
    else
      return 1;
  }
  
  double Input::potentialBeta() const {
    if (data["Physics"]["Potential"]["Beta"])
      return data["Physics"]["Potential"]["Beta"].as<double>();
    else
      return 1;
  }
  
	// ### Physics/Model ###
	
  string Input::dynamicsName() const {return data["Physics"]["Model"]["Name"].as<string>();}
  
  double Input::gamma() const {return data["Physics"]["Model"]["Gamma"].as<double>();}
  
  double Input::temperature(size_t /*iOfReplica*/) const 
  {
    if (data["Physics"]["Model"]["Temperature"])
      return data["Physics"]["Model"]["Temperature"].as<double>();
    else if (data["Physics"]["Model"]["Beta"])
      return 1 / data["Physics"]["Model"]["Beta"].as<double>();
    else if (data["Physics"]["Model"]["TemperatureLeft"] && data["Physics"]["Model"]["TemperatureRight"])
      return (temperatureLeft() + temperatureRight()) / 2;
    else if (data["Physics"]["Model"]["BetaLeft"] && data["Physics"]["Model"]["BetaRight"])
      return .5/betaLeft() + .5/betaRight();
    else 
      {cout << "Temperature not precised !" << endl;exit(1);}
  }
  
  double Input::temperatureLeft(size_t /*iOfReplica*/) const 
  {
    return data["Physics"]["Model"]["TemperatureLeft"].as<double>();
  }
  
    double Input::temperatureRight(size_t /*iOfReplica*/) const 
  {
    return data["Physics"]["Model"]["TemperatureRight"].as<double>();
  }
  
  double Input::beta(size_t /*iOfReplica*/) const 
  {
    if (data["Physics"]["Model"]["Beta"])
      return data["Physics"]["Model"]["Beta"].as<double>();
    else if (data["Physics"]["Model"]["Temperature"])
      return 1 / data["Physics"]["Model"]["Temperature"].as<double>();
    else if (data["Physics"]["Model"]["BetaLeft"] && data["Physics"]["Model"]["BetaRight"])
      return 2. / (1/betaLeft() + 1/betaRight());
    else if (data["Physics"]["Model"]["TemperatureLeft"] && data["Physics"]["Model"]["TemperatureRight"])
      return 2/(temperatureLeft() + temperatureRight());
    else
      {cout << "Beta not precised !" << endl;exit(1);}
  }
  
  double Input::betaLeft(size_t /*iOfReplica*/) const 
  {
    if (data["Physics"]["Model"]["BetaLeft"])
      return data["Physics"]["Model"]["BetaLeft"].as<double>();
    else if (data["Physics"]["Model"]["TemperatureLeft"])
      return 1 / data["Physics"]["Model"]["TemperatureLeft"].as<double>();
    else assert(false);
  }
  
  double Input::betaRight(size_t /*iOfReplica*/) const 
  {
    if (data["Physics"]["Model"]["BetaRight"])
      return data["Physics"]["Model"]["BetaRight"].as<double>();
    else if (data["Physics"]["Model"]["TemperatureRight"])
      return 1 / data["Physics"]["Model"]["TemperatureRight"].as<double>();
    else assert(false);
  }
  
  double Input::externalForceMin() const {return data["Physics"]["Model"]["Force"][0].as<double>();}
  double Input::externalForceMax() const {return data["Physics"]["Model"]["Force"][1].as<double>();}
  
  double Input::externalForce(size_t iOfReplica) const {
    if (data["Physics"]["Model"]["Force"])
      if (data["Physics"]["Model"]["Force"].size() == 2)
	return externalForceMin() + iOfReplica * (externalForceMax() - externalForceMin()) / max(1., nbOfReplicas()-1.);
      else
	return data["Physics"]["Model"]["Force"].as<double>();
    else return 0;   
  }
  
  double Input::tauBending() const
  {
    if (data["Physics"]["Model"]["Tau"])
      return data["Physics"]["Model"]["Tau"].as<double>();
    else
      return 0;
  }

  
  double Input::xi() const
  {
    if (data["Physics"]["Model"]["Xi"])
      return data["Physics"]["Model"]["Xi"].as<double>();
    else
      return 0;
  }
  
  int Input::seed() const
  {
		if (data["Physics"]["Model"]["Seed"])
			return data["Physics"]["Model"]["Seed"].as<int>();
		else
			return 0;
	}
	
	double Input::eta() const
	{
		if (dynamicsName() == "BoundaryLangevin")
			return (temperatureLeft() - temperatureRight())/2;
		else
			return externalForce();
	}
  
  
	
	//### Output ###
  
  size_t Input::decorrelationNbOfIterations(size_t iOfReplica) const
  {
    if (data["Output"]["DecorrelationTime"])
      return data["Output"]["DecorrelationTime"].as<double>() / timeStep(iOfReplica);
    else 
			return decorrelationTime(iOfReplica) / timeStep(iOfReplica);
  }
  
  double Input::decorrelationTime(size_t /*iOfReplica*/) const
  {
    if (data["Output"]["DecorrelationTime"])
      return data["Output"]["DecorrelationTime"].as<double>();
    else return 3*nbOfParticles();
  }
  
  size_t Input::outputPeriodNbOfIterations(size_t iOfReplica) const {
    if (data["Output"]["Period"])
      return data["Output"]["Period"].as<double>() / timeStep(iOfReplica);
    else
      return nbOfIterations() / 1000;
  }
  
  double Input::outputPeriodTime(size_t iOfReplica) const {
    if (data["Output"]["Period"])
      return data["Output"]["Period"].as<double>();
    else
      return outputPeriodNbOfIterations(iOfReplica) * timeStep(iOfReplica);
  }
  
  size_t Input::outputProfilePeriodNbOfIterations(size_t iOfReplica) const {
    if (data["Output"]["ProfilePeriod"])
      return data["Output"]["ProfilePeriod"].as<double>() / timeStep(iOfReplica);
    else
      return nbOfIterations() / 100;
  }
  
  double Input::outputProfilePeriodTime(size_t iOfReplica) const {
    if (data["Output"]["ProfilePeriod"])
      return data["Output"]["ProfilePeriod"].as<double>();
    else
      return outputProfilePeriodNbOfIterations(iOfReplica) * timeStep(iOfReplica);
  }
  
  /// Contains the number of values in an autocorrelation ProfilePeriod
  /// /!\ Causes a memory crash for the larger chains if too big
  int Input::nbOfAutocoPts() const
  {
		return min(1000, (int)decorrelationNbOfIterations());
	}
	
	bool Input::doFinalFlow() const
	{
		if (data["Output"]["doFinalFlow"])
			if (data["Output"]["doFinalFlow"].as<string>() == "No")
				return false;
		return true;
	}
	
	bool Input::doFinalVelocity() const
	{
		if (data["Output"]["doFinalVelocity"])
			if (data["Output"]["doFinalVelocity"].as<string>() == "No")
				return false;
		return true;
	}
	
	
	//### ControlVariate ###
  
  string Input::controlVariateName() const
  {
    if (data["ControlVariate"]["Name"])
      return data["ControlVariate"]["Name"].as<string>();
    else 
      return "None";
  }
  
  string Input::controlVariateCoeffsPath() const
  {
		if (data["ControlVariate"]["CoeffsPath"])
      return data["ControlVariate"]["CoeffsPath"].as<string>();
    else 
			assert(false);
      return "None";
  }
  
  //### Galerkin ###
  
	bool Input::doGalerkinCV() const
	{
		if (data["ControlVariate"])
			if (data["Galerkin"]["Basis"] && !data["ControlVariate"]["CoeffsPath"])
				return true;
		return false;
	}
  
  bool Input::isGalerkin() const
  {
		if (data["Galerkin"])
			if (data["Galerkin"]["Resolve"])
				if (data["Galerkin"]["Resolve"].as<string>() == "Yes")
					return true;

		return false;
	}
	
	string Input::galerkinElts() const
	{
		if (data["Galerkin"])
			if (data["Galerkin"]["Basis"])
				if (data["Galerkin"]["Basis"]["Elements"])
					return data["Galerkin"]["Basis"]["Elements"].as<string>();
		return "None";
	}
 
	size_t Input::nbOfFourier() const
	{
		if (data["Galerkin"])
			if (data["Galerkin"]["Basis"])
				if (data["Galerkin"]["Basis"]["Fourier"])
					return data["Galerkin"]["Basis"]["Fourier"].as<size_t>();

		return 0;
	}
	
	size_t Input::nbOfHermite() const
	{
		if(data["Galerkin"])
			if(data["Galerkin"]["Basis"])
				if(data["Galerkin"]["Basis"]["Hermite"])
					return data["Galerkin"]["Basis"]["Hermite"].as<size_t>();

		return 0;
	}
		
}