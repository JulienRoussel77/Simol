
#include <cstdlib>

#include "Particle.hpp"



#include<yaml-cpp/yaml.h>

int main(int argc, char* argv[])
{

  //=====================
  // COMMAND LINE PARSING
  //=====================

  YAML::Node parameters = YAML::LoadFile(argv[1]);

  //=========
  // GEOMETRY
  //=========
  
  YAML::Node geometry = parameters["Geometry"];
  size_t dimension = geometry["Dimension"].as<size_t>();
  size_t length = geometry["Length"].as<size_t>();

  //=====
  // MESH
  //=====

  YAML::Node mesh = parameters["Mesh"]["Time"];
  double timeStep = mesh["Step"].as<double>();
  size_t numberOfInstants = mesh["Instants"].as<size_t>();

  //========
  // PHYSICS
  //========

  YAML::Node physics = parameters["Physics"];
  double mass = physics["Particle"]["Mass"].as<double>();
  double initial_speed = physics["Particle"]["Speed"].as<double>();
  double initial_position = physics["Particle"]["Position"].as<double>();
  double potential_parameter = physics["Potential"]["Parameter"].as<double>();

  //==========
  // LIBRARIES
  //==========
  
  std::string linalg = parameters["Libraries"]["LinearAlgebra"].as<std::string>();

  //============
  // COMPUTATION
  //============
  
  simol::Vector<double> positions(numberOfInstants);
  simol::Vector<double> speeds(numberOfInstants);

  positions[0] = initial_position;
  speeds[0] = initial_speed;

  Particle<double> particle(mass, positions, speeds);
  simol::Potential<double> potential(potential_parameter, 2*M_PI/length);

  verlet_scheme(particle,potential,timeStep,numberOfInstants);

  std::ofstream positionFile("positions.txt");
  positionFile << particle.positions();

  std::ofstream speedFile("speeds.txt");
  speedFile << particle.speeds();
  
  return EXIT_SUCCESS;
}
