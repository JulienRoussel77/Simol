
#include <cstdlib>

#include "ParticleSystem.hpp"
#include "HamiltonDynamics.hpp"


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
  YAML::Node parameters = YAML::LoadFile(cmd.inputFileName());

  YAML::Node geometry = parameters["Geometry"];
  size_t dimension = geometry["Dimension"].as<size_t>();
  size_t length = geometry["Length"].as<size_t>();

  YAML::Node mesh = parameters["Mesh"]["Time"];
  double timeStep = mesh["Step"].as<double>();
  size_t numberOfInstants = mesh["Number"].as<size_t>();

  YAML::Node physics = parameters["Physics"];
  double mass = physics["Particle"]["Mass"].as<double>();
  double initial_speed = physics["Particle"]["Speed"].as<double>();
  double initial_position = physics["Particle"]["Position"].as<double>();
  double potential_parameter = physics["Potential"]["Parameter"].as<double>();

  std::string outputFilename = parameters["Output"]["Filename"].as<std::string>();

  //============
  // COMPUTATION
  //============

  size_t numberOfParticles = 1;

  simol::ParticleSystem<double> system(numberOfParticles, mass, initial_position, initial_speed);
  simol::Potential<double> potential(potential_parameter, 2*M_PI/length);

  simol::HamiltonDynamics<double> model(mass, potential);
  
  std::ofstream outputFile(outputFilename);
  
  for (size_t instantIndex  =1; instantIndex < numberOfInstants; ++instantIndex)
  {
    double instant = instantIndex * timeStep;
    system.simulate(timeStep, model, outputFile);
  }

  return EXIT_SUCCESS;
}
