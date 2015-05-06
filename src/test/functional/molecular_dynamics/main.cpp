
#include <cstdlib>

#include "ParticleSystem.hpp"
#include "HamiltonianDynamics.hpp"


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
  double finalInstant = mesh["Final"].as<double>();

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
  size_t numberOfInstants = std::floor(finalInstant/timeStep)+2;

  simol::ParticleSystem<double> system(numberOfParticles, mass, initial_position, initial_speed);
  simol::Potential<double> potential(potential_parameter, 2*M_PI/length);

  simol::HamiltonianDynamics<double> model(mass, potential);
  
  std::ofstream outputFile(outputFilename);
  
  for (double instant = timeStep; instant < finalInstant; instant+=timeStep)
  {
    system.simulate(instant, model, outputFile);
   /* for (auto&& particle : system.configuration())
    {
      verlet(particle,potential,timeStep);
      outputFile << instant << " " << particle.position() << " " << particle.speed() << std::endl;
    }*/
  }

  return EXIT_SUCCESS;
}
