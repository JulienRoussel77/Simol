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

  std::string Input::systemName() const {return data["Physics"]["System"]["Name"].as<std::string>();}
  size_t Input::numberOfParticles() const {return 1;}
  double Input::mass() const {return data["Physics"]["System"]["Mass"].as<double>();} 	
  double Input::initialPosition(int const& i) const {return data["Physics"]["System"]["Position"].as<double>();}       
  double Input::initialSpeed(int const& i) const {return data["Physics"]["System"]["Speed"].as<double>();}
  
  int Input::numberOfReplica() const {return data["Physics"]["Replica"]["Number"].as<int>();}

  std::string Input::outputFilename() const {return data["Output"]["Filename"].as<std::string>();}

}