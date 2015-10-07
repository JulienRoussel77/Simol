#ifndef SIMOL_INPUT_HPP
#define SIMOL_INPUT_HPP

#include<yaml-cpp/yaml.h>
#include <getopt.h>
#include "Vector.hpp"

#include "io/CommandLine.hpp"

#include <string>
#include <iostream>

namespace simol{
  
  class Input{
    YAML::Node data;
  public:
    Input(CommandLine cmd);
    int dimension() const;
    double length() const;
    
    double timeStep() const;
    size_t numberOfIterations() const;
    
    std::string potentialName() const;
    
    //Sinus
    double amplitude() const;
    
    //DoubleWell
    double height() const;
    double interWell() const;
    
    //Harmonic
    double stiffness() const;
	
    std::string dynamicsName() const;
    double gamma() const;   
    double temperature() const;
    double force() const;
    
    std::string systemName() const;
    size_t numberOfParticles() const;
    double mass() const;
    

      
    double initialPosition(int const& i) const;
    
    double initialSpeed(int const& i) const;
    
    int numberOfReplica() const;
    
    std::string outputFoldername() const;

  };
  

}

#endif