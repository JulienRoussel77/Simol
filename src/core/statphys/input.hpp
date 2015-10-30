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
    
    double timeStepMin() const;
    double timeStepMax() const;
    double timeStep(size_t indexOfReplica=0) const;
    size_t numberOfIterations(size_t indexOfReplica=0) const;
    
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
    double beta() const;
    double betaLeft() const;
    double betaRight() const;
    //bool externalForceVarying() const;
    double externalForceMin() const;
    double externalForceMax() const;
    double externalForce(size_t indexOfReplica=0) const;
    
    std::string systemName() const;
    size_t numberOfParticles() const;
    double mass() const;
    

      
    double initialPosition(int const& i=0) const;
    
    double initialMomentum(int const& i=0) const;
    
    size_t numberOfReplicas() const;
    
    //bool doComputeCorrelations() const;
    size_t decorrelationNumberOfIterations(size_t indexOfReplica=0) const;
    double decorrelationTime(size_t indexOfReplica=0) const;
    std::string outputFoldername() const;
    size_t outputPeriodNumberOfIterations(size_t indexOfReplica=0) const;
    double outputPeriodTime(size_t indexOfReplica=0) const;
    
    std::string controlVariateName() const;
  };
  

}

#endif