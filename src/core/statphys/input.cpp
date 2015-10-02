#include "input.hpp"

namespace simol {
  

  Input::Input(CommandLine cmd):data(YAML::LoadFile(cmd.inputFileName()))
  {
    std::cout << "Input read in " << cmd.inputFileName() << std::endl;
  }

  int Input::dimension() const {return data["Geometry"]["Dimension"].as<int>();}

  double Input::length() const {return data["Geometry"]["Length"].as<double>();}

  double Input::timeStep() const {return data["Mesh"]["Time"]["Step"].as<double>();}

  size_t Input::numberOfIterations() const {return data["Mesh"]["Time"]["Number"].as<double>();}

  double Input::potParameter() const {return data["Physics"]["Potential"]["Parameter"].as<double>();}

  std::string Input::name() const {return data["Physics"]["Model"]["Name"].as<std::string>();}

  double Input::gamma() const {return data["Physics"]["Model"]["Gamma"].as<double>();}

  double Input::temperature() const {return data["Physics"]["Model"]["Temperature"].as<double>();}

  size_t Input::numberOfParticles() const {return 1;}

  double Input::mass() const {return data["Physics"]["Particle"]["Mass"].as<double>();} 
	
  double Input::initialPosition(int const& i) const {return data["Physics"]["Particle"]["Position"].as<double>();} 
      
  double Input::initialSpeed(int const& i) const {return data["Physics"]["Particle"]["Speed"].as<double>();}

  std::string Input::outputFilename() const {return data["Output"]["Filename"].as<std::string>();}

}