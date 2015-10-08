#include "input.hpp"

using std::cout; 
using std::endl; 

namespace simol {
  

  Input::Input(CommandLine cmd):data(YAML::LoadFile(cmd.inputFileName()))
  {
    std::cout << "Input read in " << cmd.inputFileName() << std::endl;
  }

  int Input::dimension() const {return data["Geometry"]["Dimension"].as<int>();}

  double Input::length() const {return data["Geometry"]["Length"].as<double>();}

  double Input::timeStep() const {return data["Mesh"]["Time"]["Step"].as<double>();}

  size_t Input::numberOfIterations() const {return data["Mesh"]["Time"]["Number"].as<double>();}
  
  std::string Input::potentialName() const {return data["Physics"]["Potential"]["Name"].as<std::string>();}

  //Sinusoidal
  double Input::amplitude() const {return data["Physics"]["Potential"]["Amplitude"].as<double>();}
  
  //DoubleWell
  double Input::height() const {return data["Physics"]["Potential"]["Height"].as<double>();}
  double Input::interWell() const {return data["Physics"]["Potential"]["InterWell"].as<double>();}
  
  //Harmonic
  double Input::stiffness() const {return data["Physics"]["Potential"]["Stiffness"].as<double>();}

  std::string Input::dynamicsName() const {return data["Physics"]["Model"]["Name"].as<std::string>();}
  double Input::gamma() const {return data["Physics"]["Model"]["Gamma"].as<double>();}
  double Input::temperature() const {return data["Physics"]["Model"]["Temperature"].as<double>();}
  

  
  bool Input::externalForceVarying() const {
    if (data["Physics"]["Model"]["ForceMin"] && data["Physics"]["Model"]["ForceMax"] && !data["Physics"]["Model"]["Force"]) return true;
    else if (!data["Physics"]["Model"]["ForceMin"] && !data["Physics"]["Model"]["ForceMax"] && data["Physics"]["Model"]["Force"]) return false;
    else {cout << "External force input incoherent !" << endl;exit(1);}
  }
  double Input::externalForceMin() const {return data["Physics"]["Model"]["ForceMin"].as<double>();}
  double Input::externalForceMax() const {return data["Physics"]["Model"]["ForceMax"].as<double>();}
  
  double Input::externalForce(int const& indexOfReplica) const {
    if (!externalForceVarying())
      if (data["Physics"]["Model"]["Force"])
	return data["Physics"]["Model"]["Force"].as<double>();
      else return 0;
    else
      return externalForceMin() + indexOfReplica * (externalForceMax() - externalForceMin()) / numberOfReplicas();
  }

  std::string Input::systemName() const {return data["Physics"]["System"]["Name"].as<std::string>();}
  size_t Input::numberOfParticles() const {return 1;}
  
  double Input::mass() const {
    if (data["Physics"]["System"]["Mass"])
      return data["Physics"]["System"]["Mass"].as<double>();
    else return 1;
  } 	
  
  double Input::initialPosition(int const& i) const {
    if (data["Physics"]["System"]["Position"])
      return data["Physics"]["System"]["Position"].as<double>();
    else return 0;
  }   
  
    double Input::initialMomentum(int const& i) const {
    if (data["Physics"]["System"]["Momentum"])
      return data["Physics"]["System"]["Momentum"].as<double>();
    else return 0;
  }  
  
  size_t Input::numberOfReplicas() const {
    if (data["Physics"]["Replica"]["Number"])
      return data["Physics"]["Replica"]["Number"].as<int>();
    else return 1;
  }

  //std::string Input::outputFilename() const {return data["Output"]["Filename"].as<std::string>();}
  std::string Input::outputFoldername() const {return "../output/"+dynamicsName()+"/"+systemName()+"/"+potentialName()+"/";}

}