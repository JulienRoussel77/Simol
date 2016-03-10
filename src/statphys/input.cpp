#include "input.hpp"

namespace simol {


  Input::Input(CommandLine cmd):
		data(YAML::LoadFile(cmd.inputFileName())),
		inputPath_(cmd.inputFileName()),
		inputFlux_(inputPath())
		//inputSettings_(settingsPath()),
		//settingsSize_(settingsSize())
  {
		/*std::cout << "SettingsPath = " << settingsPath() << std::endl;
		assert(inputSettings_.is_open());
        std::cout << "SettingsSize = " << settingsSize_ << std::endl;*/

    if (data["Physics"]["System"]["Position"])
      if (data["Physics"]["System"]["Position"].size() == 2)
      {
				positionMin_ = data["Physics"]["System"]["Position"][0].as<double>();
				positionMax_ = data["Physics"]["System"]["Position"][1].as<double>();
      }


    /*if (settingsPath() != "" && settingsSize_ != -1)
		{
			settingsPositions_ = vector<double>(settingsSize_);
			settingsMomenta_ = vector<double>(settingsSize_);
			for (size_t iOfParticle=0; (int)iOfParticle < settingsSize_; iOfParticle++)
			{
				readItem(inputSettings_);
				assert( (int) readItem(inputSettings_) == (int) iOfParticle);
				settingsPositions_[iOfParticle] = readItem(inputSettings_);
				settingsMomenta_[iOfParticle] = readItem(inputSettings_);
				for (int i=0; i<4; i++)
					readItem(inputSettings_);
                    std::cout << iOfParticle << " " << settingsPositions_[iOfParticle] << " " << settingsMomenta_[iOfParticle] << std::endl;
			}

			//Estimating the tauBending
			//double estimTau = 2 * (settingsPositions_[settingsSize()-1] - 2*settingsPositions_[(settingsSize()-1)/2] + settingsPositions_[(settingsSize()-1)/2]
			double estimTau = settingsPositions_[settingsSize()-1] / pow(settingsSize()-1, 4);
            std::cout << "Estimated bending = " << estimTau << std::endl;

            std::cout << "file read" << std::endl;

			initialPositions_ = vector<double>(numberOfParticles());
			initialMomenta_ = vector<double>(numberOfParticles());
			double ratio = settingsSize_ / (double) numberOfParticles();
			for (size_t iOfParticle=0; iOfParticle < numberOfParticles(); iOfParticle++)
			{
				double fracIndex = ratio * iOfParticle;
				int iInf = (int) (fracIndex);
				int iSup = ceil  (fracIndex);
				double alpha = fracIndex - iInf;
				assert(0 <= alpha && alpha <= 1);
				initialPositions_[iOfParticle] = pow(ratio, -4) *
						((1-alpha) * settingsPositions_[iInf] + alpha * settingsPositions_[iSup]);
				initialMomenta_[iOfParticle]   = (1-alpha) * settingsMomenta_[iInf]   + alpha * settingsMomenta_[iSup];

                std::cout << iOfParticle << " " << alpha << " " << settingsPositions_[iInf] << " " << settingsPositions_[iSup] << " " << initialPositions_[iOfParticle] << std::endl;
			}
		}*/
  }

  const std::string& Input::inputPath() const
  {
		return inputPath_;
	}

	const std::ifstream& Input::inputFlux() const
	{
		return inputFlux_;
	}

  int Input::dimension() const {return data["Geometry"]["Dimension"].as<int>();}

  double Input::length() const {return data["Geometry"]["Length"].as<double>();}

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

  double Input::timeStep(size_t indexOfReplica) const
  {
    //return timeStepMin() + indexOfReplica * (timeStepMax() - timeStepMin()) / numberOfReplicas();

    return timeStepMin() * pow(timeStepMax() / timeStepMin(), indexOfReplica / std::max(1., numberOfReplicas()-1.));
  }

  //double Input::timeStep() const {return data["Mesh"]["Time"]["Step"].as<double>();}

  size_t Input::numberOfIterations(size_t indexOfReplica) const
  {
    if (data["Mesh"]["Time"]["Number"])
      return data["Mesh"]["Time"]["Number"].as<double>();
    else if (data["Mesh"]["Time"]["FinalTime"])
      return data["Mesh"]["Time"]["FinalTime"].as<double>() / timeStep(indexOfReplica);
    else
    {std::cout << "Number of Iterations not specified !" << std::endl;exit(1);}
  }

  std::string Input::potentialName() const {return data["Physics"]["Potential"]["Name"].as<std::string>();}

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


  std::string Input::dynamicsName() const {return data["Physics"]["Model"]["Name"].as<std::string>();}
  double Input::gamma() const {return data["Physics"]["Model"]["Gamma"].as<double>();}

  double Input::temperature(size_t /*indexOfReplica*/) const
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
      {std::cout << "Temperature not precised !" << std::endl;exit(1);}
  }

  double Input::temperatureLeft(size_t /*indexOfReplica*/) const
  {
    return data["Physics"]["Model"]["TemperatureLeft"].as<double>();
  }

    double Input::temperatureRight(size_t /*indexOfReplica*/) const
  {
    return data["Physics"]["Model"]["TemperatureRight"].as<double>();
  }

  double Input::beta(size_t /*indexOfReplica*/) const
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
      {std::cout << "Beta not precised !" << std::endl;exit(1);}
  }

  double Input::betaLeft(size_t /*indexOfReplica*/) const
  {
    if (data["Physics"]["Model"]["BetaLeft"])
      return data["Physics"]["Model"]["BetaLeft"].as<double>();
    else if (data["Physics"]["Model"]["TemperatureLeft"])
      return 1 / data["Physics"]["Model"]["TemperatureLeft"].as<double>();
    else assert(false);
  }

  double Input::betaRight(size_t /*indexOfReplica*/) const
  {
    if (data["Physics"]["Model"]["BetaRight"])
      return data["Physics"]["Model"]["BetaRight"].as<double>();
    else if (data["Physics"]["Model"]["TemperatureRight"])
      return 1 / data["Physics"]["Model"]["TemperatureRight"].as<double>();
    else assert(false);
  }

  /*bool Input::externalForceVarying() const {
    if (data["Physics"]["Model"]["ForceMin"] && data["Physics"]["Model"]["ForceMax"] && !data["Physics"]["Model"]["Force"]) return true;
    else if (!data["Physics"]["Model"]["ForceMin"] && !data["Physics"]["Model"]["ForceMax"]) return false;
    else {std::cout << "External force input incoherent !" << std::endl;exit(1);}
  }*/

  double Input::externalForceMin() const {return data["Physics"]["Model"]["Force"][0].as<double>();}
  double Input::externalForceMax() const {return data["Physics"]["Model"]["Force"][1].as<double>();}

  double Input::externalForce(size_t indexOfReplica) const {
    if (data["Physics"]["Model"]["Force"])
      if (data["Physics"]["Model"]["Force"].size() == 2)
	return externalForceMin() + indexOfReplica * (externalForceMax() - externalForceMin()) / std::max(1., numberOfReplicas()-1.);
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

  std::string Input::systemName() const {return data["Physics"]["System"]["Name"].as<std::string>();}

  size_t Input::numberOfParticles() const {
    if (data["Physics"]["System"]["Number"])
      return data["Physics"]["System"]["Number"].as<size_t>();
    else return 1;
  }

  double Input::mass() const {
    if (data["Physics"]["System"]["Mass"])
      return data["Physics"]["System"]["Mass"].as<double>();
    else return 1;
  }

  double Input::initialPosition(int const& i) const {
    if (data["Physics"]["System"]["Position"])
		{
      if (data["Physics"]["System"]["Position"].size() == 1)
				return data["Physics"]["System"]["Position"].as<double>();
      else
				return positionMin_ + (i + .5)/numberOfParticles() * (positionMax_ - positionMin_);
		}
		/*else if (settingsSize_ != -1)
		{
        std::cout << "using settings for q : " << i << "->" << initialPositions_[i] << std::endl;
			return initialPositions_[i];
		}*/

		else return 0;
  }

  double Input::initialMomentum(int const& /*i*/) const
  {
    if (data["Physics"]["System"]["Momentum"])
      return data["Physics"]["System"]["Momentum"].as<double>();
		/*else if (settingsSize_ != -1)
		{
        std::cout << "using settings for p : " << i << "->" << initialMomenta_[i] << std::endl;
			return initialMomenta_[i];
		}*/
    else return 0;
    //else return (i < numberOfParticles()/2)?.2:-.2;
  }

  std::string Input::settingsPath() const
  {
		if (data["Physics"]["System"]["SettingsPath"])
			return finalOutputFoldername()+data["Physics"]["System"]["SettingsPath"].as<std::string>();
		else return "";
	}

	int Input::settingsSize() const
  {
		if (data["Physics"]["System"]["SettingsSize"])
			return data["Physics"]["System"]["SettingsSize"].as<int>();
		else return -1;
	}


  size_t Input::numberOfReplicas() const {
    if (data["Physics"]["Replicas"]["Number"])
      return data["Physics"]["Replicas"]["Number"].as<int>();
    else return 1;
  }

  /*bool Input::doComputeCorrelations() const {
    return data["Output"]["DecorrelationTime"];
  }*/

  size_t Input::decorrelationNumberOfIterations(size_t indexOfReplica) const
  {
    if (data["Output"]["DecorrelationTime"])
      return data["Output"]["DecorrelationTime"].as<double>() / timeStep(indexOfReplica);
    else return 0;
  }

  double Input::decorrelationTime(size_t /*indexOfReplica*/) const
  {
    if (data["Output"]["DecorrelationTime"])
      return data["Output"]["DecorrelationTime"].as<double>();
    else return 0;
  }

  //std::string Input::outputFilename() const {return data["Output"]["Filename"].as<std::string>();}
  std::string Input::outputFoldername() const {
    //std::cout << "../output/"+dynamicsName()+"/"+systemName()+"/"+potentialName()+"/" << std::endl;
      std::string foldername = "../output/"+dynamicsName()+"/"+systemName()+"/"+potentialName()+"/";
    if (controlVariateName() != "None")
      foldername += controlVariateName()+"/";
    if (data["Output"]["Foldername"])
      foldername += data["Output"]["Foldername"].as<std::string>()+"/";
    return foldername;
  }

  std::string Input::finalOutputFoldername() const {
    //std::cout << "../output/"+dynamicsName()+"/"+systemName()+"/"+potentialName()+"/" << std::endl;
      std::string finalFoldername = "../output/"+dynamicsName()+"/"+systemName()+"/"+potentialName()+"/";
    if (controlVariateName() != "None")
      finalFoldername += controlVariateName()+"/";
    return finalFoldername;
  }

  size_t Input::outputPeriodNumberOfIterations(size_t indexOfReplica) const {
    if (data["Output"]["Period"])
      return data["Output"]["Period"].as<double>() / timeStep(indexOfReplica);
    else
      return 1;
  }

  double Input::outputPeriodTime(size_t /*indexOfReplica*/) const {
    if (data["Output"]["Period"])
      return data["Output"]["Period"].as<double>();
    else
      return 1;
  }

  std::string Input::controlVariateName() const
  {
    if (data["ControlVariate"]["Name"])
      return data["ControlVariate"]["Name"].as<std::string>();
    else
      return "None";
  }

  //####Galerkin####

  	size_t Input::numberOfFourier() const
  	{
			if(data["Galerkin"])
				if(data["Galerkin"]["Basis"])
					if(data["Galerkin"]["Basis"]["Fourier"])
						return data["Galerkin"]["Basis"]["Fourier"].as<size_t>();

			return 0;
		}

		size_t Input::numberOfHermite() const
		{
			if(data["Galerkin"])
				if(data["Galerkin"]["Basis"])
					if(data["Galerkin"]["Basis"]["Hermite"])
						return data["Galerkin"]["Basis"]["Hermite"].as<size_t>();

			return 0;
		}



		double Input::readItem(std::ifstream& in)
		{
			double item;
			in >> item;
			return item;
		}

}
