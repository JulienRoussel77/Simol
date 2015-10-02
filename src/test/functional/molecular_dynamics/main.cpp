
#include <cstdlib>
#include <iostream>

#include "particleSystem.hpp"
#include "dynamics.hpp"
#include "Vector.hpp"
#include "input.hpp"

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
    
  simol::Input input(cmd);

  //============
  // COMPUTATION
  //============

  /*size_t numberOfParticles = 1;
  
  std::vector<double> masses(numberOfParticles);
  for (int i=0; i<numberOfParticles; i++)
    masses[i] = mass;
  
  std::vector<dvec> positions(numberOfParticles);
  for (int i=0; i<numberOfParticles; i++)
    positions[i] = initial_position;
  
  std::vector<dvec> momenta(numberOfParticles);
  for (int i=0; i<numberOfParticles; i++)
    momenta[i] = initial_speed;

  simol::ParticleSystem system(numberOfParticles, masses, positions, momenta);*/
  
  simol::ParticleSystem system(input);
  
  simol::Potential potential(input);

  simol::Dynamics* model = simol::createDynamics(input);
  
  std::ofstream outputFile(input.outputFilename());
  
  std::cout << "Output written in " << input.outputFilename() << std::endl;
  
  for (size_t instantIndex  =1; instantIndex < input.numberOfIterations(); ++instantIndex)
  {
    double instant = instantIndex * input.timeStep();
    system.simulate(input.timeStep(), model);
  }

  //===========
  // ATOM CHAIN
  //===========

  std::cout << "Fin de la simulation" << std::endl;

  return EXIT_SUCCESS;
}
