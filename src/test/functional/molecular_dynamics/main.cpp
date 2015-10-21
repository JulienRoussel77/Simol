
#include <cstdlib>
#include <iostream>

#include "multiSystem.hpp"
#include "particleSystem.hpp"
#include "dynamics.hpp"
#include "Vector.hpp"
#include "input.hpp"

#include <cmath>

#include<yaml-cpp/yaml.h>
#include <getopt.h>

#include "io/CommandLine.hpp"

using std::cout; 
using std::endl; 

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
  
  simol::MultiSystem replica(input);

  //simol::Dynamics* model = simol::createDynamics(input);
  
  replica.launch();
  

  

  std::cout << "Fin de la simulation" << std::endl;

  return EXIT_SUCCESS;
}
