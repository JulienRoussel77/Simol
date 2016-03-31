#include <cstdlib>
#include <iostream>

#include "simulation.hpp"
#include "system.hpp"
#include "Dynamics.hpp"
#include "Vector.hpp"
#include "input.hpp"
#include "tools.hpp"

#include <cmath>

#include<yaml-cpp/yaml.h>
#include <getopt.h>

#include "CommandLine.hpp"
#include <time.h>

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
  double totalTime = clock();

  //=====================
  // COMMAND LINE PARSING
  //=====================

  simol::CommandLine cmd(argc,argv);

  //===================
  // INPUT FILE LOADING
  //===================

  cout << "Input read in " << cmd.inputFileName() << endl;
  simol::Input input(cmd);

  //============
  // COMPUTATION
  //============

  if (input.isGalerkin())
    throw std::invalid_argument("The input must correspond to a MD simulation !");

  simol::Simulation<simol::DPDE, simol::Isolated> simu(input);
  simu.launch();

  displayTime(clock() - totalTime);

  return EXIT_SUCCESS;
}
