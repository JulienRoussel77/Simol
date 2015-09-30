
#include <cstdlib>
#include <iostream>

#include "particleSystem.hpp"
#include "dynamics.hpp"
#include "atomChain.hpp"

#include <cmath>

#include<yaml-cpp/yaml.h>
#include <getopt.h>

#include "io/CommandLine.hpp"

int main(int argc, char* argv[])
{

  //=====================
  // COMMAND LINE PARSING
  //=====================

  simol::CommandLine cmd(argc,argv);

  //===================
  // INPUT FILE LOADING
  //===================
  
  std::cout << "Input read in " << cmd.inputFileName() << std::endl;
    
  YAML::Node input = YAML::LoadFile(cmd.inputFileName());

  YAML::Node geometry = input["Geometry"];
  size_t dimension = geometry["Dimension"].as<size_t>();
  size_t length = geometry["Length"].as<size_t>();

  YAML::Node mesh = input["Mesh"]["Time"];
  double timeStep = mesh["Step"].as<double>();
  size_t numberOfInstants = mesh["Number"].as<size_t>();

  YAML::Node physics = input["Physics"];
  double mass = physics["Particle"]["Mass"].as<double>();
  double initial_speed = physics["Particle"]["Speed"].as<double>();
  double initial_position = physics["Particle"]["Position"].as<double>();
  double potential_parameter = physics["Potential"]["Parameter"].as<double>();

  std::string outputFilename = input["Output"]["Filename"].as<std::string>();

  //============
  // COMPUTATION
  //============

  size_t numberOfParticles = 1;

  simol::ParticleSystem system(numberOfParticles, mass, initial_position, initial_speed);
  simol::Potential potential(potential_parameter, 2*M_PI/length);

  simol::Dynamics* model = simol::createDynamics(input);
  
  std::ofstream outputFile(outputFilename);
  
  std::cout << "Output written in " << outputFilename << std::endl;
  
  for (size_t instantIndex  =1; instantIndex < numberOfInstants; ++instantIndex)
  {
    double instant = instantIndex * timeStep;
    system.simulate(timeStep, model, outputFile);
  }

  //===========
  // ATOM CHAIN
  //===========
  
  size_t numberOfAtoms = 3;
  double leftPosition = 0;
  double leftMomentum = 0;
  double leftTemperature = 0;
  double rightTemperature = 1;
  simol::AtomChain chain(numberOfAtoms, mass, leftPosition, leftMomentum, leftTemperature, rightTemperature, potential);

  std::cout << "Fin de la simulation" << std::endl;

  return EXIT_SUCCESS;
}
