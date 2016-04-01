#include "Input.hpp"

using std::cout;
using std::endl;
using std::string;
using std::max;
using std::min;

const double defaultMass = 1;
const double defaultLength = 2 * M_PI;
const int defaultNbOfParticles = 1;
const double defaultPotentialCoeff = 1;
const double defaultSeed = 0;
const int defaultOutputPeriod = 1;         //10
const int defaultOutputProfilePeriod = 1;  //100
//-- default values for dynamics --
const double defaultHeatCapacity = 1;
const double defaultGamma = 1; 
const double defaultXi = 0;
const double defaultExternalForce = 0;
const double defaultTauBending = 0;
//-- default values for output --
const int maxNbOfAutocoPts = 1000;

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

		if (doFileSettings())
		{
			assert(inputFlux_.is_open());
			cout << "Reading the settings from " << settingsPath() << "...";
			initialPositions_ = vector<Vector<double>>(nbOfParticles());
			initialMomenta_ = vector<Vector<double>>(nbOfParticles());
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
    return "../../../output/"+dynamicsName()+"/"+systemName()+"/"+potentialName()+"/";
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
  
  /*const std::shared_ptr<RNG>& Input::rng() const
  {
    return rng_;
  }
  
  std::shared_ptr<RNG>& Input::rng()
  {
    return rng_;
  }*/

	//### Geometry ###

  int Input::dimension() const {return data["Geometry"]["Dimension"].as<int>();}

  double Input::length() const
  {
		if (data["Geometry"]["Length"])
			return data["Geometry"]["Length"].as<double>();
		else
			return defaultLength;
	}

	// ### Mesh/Time ###

  double Input::timeStep() const
  {
    if (data["Mesh"]["Time"]["Step"])
      return data["Mesh"]["Time"]["Step"].as<double>();
    throw std::invalid_argument("Timestep missing");
  }

  ///
	///Number of iterations in the "output" part of the simulation
  size_t Input::nbOfIterations() const
  {
    if (data["Mesh"]["Time"]["Number"])
      return data["Mesh"]["Time"]["Number"].as<size_t>();
    else if (data["Mesh"]["Time"]["FinalTime"])
      return data["Mesh"]["Time"]["FinalTime"].as<double>() / timeStep();
    else
			{cout << "Nb of Iterations not specified !" << endl;exit(1);}
  }

  ///
  ///Number of iterations in the thermalization part (temperatures imposed everywhere)
  size_t Input::nbOfThermalIterations() const
  {
    if (data["Mesh"]["Time"]["ThermalNumber"])
      return data["Mesh"]["Time"]["ThermalNumber"].as<size_t>();
    else if (data["Mesh"]["Time"]["ThermalTime"])
      return data["Mesh"]["Time"]["ThermalTime"].as<double>() / timeStep();
    else
			return 0;
  }

  ///
  ///Number of iterations in the burning part (no output)
  size_t Input::nbOfBurningIterations() const
  {
    if (data["Mesh"]["Time"]["BurningNumber"])
      return data["Mesh"]["Time"]["BurningNumber"].as<size_t>();
    else if (data["Mesh"]["Time"]["BurningTime"])
      return data["Mesh"]["Time"]["BurningTime"].as<double>() / timeStep();
    else
			return 0;
  }

  // ### Physics/System ###

  string Input::systemName() const {return data["Physics"]["System"]["Name"].as<string>();}

  size_t Input::nbOfParticles() const {
    if (data["Physics"]["System"]["Number"])
      return data["Physics"]["System"]["Number"].as<size_t>();
    else return defaultNbOfParticles;
  }

  double Input::mass() const {
    if (data["Physics"]["System"]["Mass"])
      return data["Physics"]["System"]["Mass"].as<double>();
    else return defaultMass;
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
  
  Vector<double> Input::initialPosition(int const& iOfParticle) const {
    Vector<double> q0(dimension(), 0);
    if (data["Physics"]["System"]["Position"])
      q0(0) = data["Physics"]["System"]["Position"].as<double>();
		else if (doFileSettings())
		{
			cout << "using settings for q : " << iOfParticle << "->" << initialPositions_[iOfParticle] << endl;
			q0 = initialPositions_[iOfParticle];
		}
    return q0;
  }   
  
  Vector<double> Input::initialMomentum(int const& iOfParticle) const 
  {
    Vector<double> p0(dimension(), 0);
    if (data["Physics"]["System"]["Momentum"])
      p0(0) =  data["Physics"]["System"]["Momentum"].as<double>();
		else if (doFileSettings())
		{
			cout << "using settings for p : " << iOfParticle << "->" << initialMomenta_[iOfParticle] << endl;
			p0 = initialMomenta_[iOfParticle];
		}
    return p0;
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
      return defaultPotentialCoeff;
  }

  //DoubleWell
  double Input::height() const {return data["Physics"]["Potential"]["Height"].as<double>();}
  double Input::interWell() const {return data["Physics"]["Potential"]["InterWell"].as<double>();}

  //Harmonic
  double Input::potentialStiffness() const {
    if (data["Physics"]["Potential"]["Stiffness"])
      return data["Physics"]["Potential"]["Stiffness"].as<double>();
    else
      return defaultPotentialCoeff;
  }

  //Quadratic
  double Input::potentialAlpha() const {
    if (data["Physics"]["Potential"]["Alpha"])
      return data["Physics"]["Potential"]["Alpha"].as<double>();
    else
      return defaultPotentialCoeff;
  }

  double Input::potentialBeta() const {
    if (data["Physics"]["Potential"]["Beta"])
      return data["Physics"]["Potential"]["Beta"].as<double>();
    else
      return defaultPotentialCoeff;
  }

	// ### Physics/Model ###

  string Input::dynamicsName() const {return data["Physics"]["Model"]["Name"].as<string>();}

  double Input::gamma() const {return data["Physics"]["Model"]["Gamma"].as<double>();}

  double Input::temperature() const
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

  double Input::temperatureLeft() const
  {
    return data["Physics"]["Model"]["TemperatureLeft"].as<double>();
  }

  double Input::temperatureRight() const
  {
    return data["Physics"]["Model"]["TemperatureRight"].as<double>();
  }

  double Input::beta() const
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
      throw std::invalid_argument("Beta not precised !");
  }

  double Input::betaLeft() const
  {
    if (data["Physics"]["Model"]["BetaLeft"])
      return data["Physics"]["Model"]["BetaLeft"].as<double>();
    else if (data["Physics"]["Model"]["TemperatureLeft"])
      return 1 / data["Physics"]["Model"]["TemperatureLeft"].as<double>();
    else assert(false);
  }

  double Input::betaRight() const
  {
    if (data["Physics"]["Model"]["BetaRight"])
      return data["Physics"]["Model"]["BetaRight"].as<double>();
    else if (data["Physics"]["Model"]["TemperatureRight"])
      return 1 / data["Physics"]["Model"]["TemperatureRight"].as<double>();
    else assert(false);
  }

  double Input::externalForce() const {
    if (data["Physics"]["Model"]["Force"])
      return data["Physics"]["Model"]["Force"].as<double>();
    else return defaultExternalForce;
  }

  double Input::tauBending() const
  {
    if (data["Physics"]["Model"]["Tau"])
      return data["Physics"]["Model"]["Tau"].as<double>();
    else
      return defaultTauBending;
  }


  double Input::xi() const
  {
    if (data["Physics"]["Model"]["Xi"])
      return data["Physics"]["Model"]["Xi"].as<double>();
    else
      return defaultXi;
  }

  int Input::seed() const
  {
    if (data["Physics"]["Model"]["Seed"])
      return data["Physics"]["Model"]["Seed"].as<int>();
    else
      return defaultSeed;
  }

  double Input::eta() const
  {
    if (dynamicsName() == "BoundaryLangevin")
      return (temperatureLeft() - temperatureRight())/2;
    else
      return externalForce();
  }
  
  double Input::heatCapacity() const 
  {
    if (data["Physics"]["Model"]["HeatCapacity"])
      return data["Physics"]["Model"]["HeatCapacity"].as<double>();
    else
      return defaultHeatCapacity;
  }
  
  //### Output ###

  size_t Input::decorrelationNbOfIterations() const
  {
    if (data["Output"]["DecorrelationTime"])
      return data["Output"]["DecorrelationTime"].as<double>() / timeStep();
    else
			return decorrelationTime() / timeStep();
  }

  double Input::decorrelationTime() const
  {
    if (data["Output"]["DecorrelationTime"])
      return data["Output"]["DecorrelationTime"].as<double>();
    else return 3*nbOfParticles();
  }

  size_t Input::outputPeriodNbOfIterations() const {
    if (data["Output"]["Period"])
      return data["Output"]["Period"].as<double>() / timeStep();
    else
      return defaultOutputPeriod;
  }

  double Input::outputPeriodTime() const {
    if (data["Output"]["Period"])
      return data["Output"]["Period"].as<double>();
    else
      return outputPeriodNbOfIterations() * timeStep();
  }

  size_t Input::outputProfilePeriodNbOfIterations() const {
    if (data["Output"]["ProfilePeriod"])
      return data["Output"]["ProfilePeriod"].as<double>() / timeStep();
    else
      return defaultOutputProfilePeriod;
  }

  double Input::outputProfilePeriodTime() const {
    if (data["Output"]["ProfilePeriod"])
      return data["Output"]["ProfilePeriod"].as<double>();
    else
      return outputProfilePeriodNbOfIterations() * timeStep();
  }

  /// Contains the number of values in an autocorrelation ProfilePeriod
  /// /!\ Causes a memory crash for the larger chains if too big
  int Input::nbOfAutocoPts() const
  {
		return min(maxNbOfAutocoPts, (int)decorrelationNbOfIterations());
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

		throw std::invalid_argument("Number of Fourier modes missing");
	}

	size_t Input::nbOfHermite() const
	{
		if(data["Galerkin"])
			if(data["Galerkin"]["Basis"])
				if(data["Galerkin"]["Basis"]["Hermite"])
					return data["Galerkin"]["Basis"]["Hermite"].as<size_t>();

    throw std::invalid_argument("Number of Hermite modes missing");
	}

}
