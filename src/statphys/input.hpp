#ifndef SIMOL_INPUT_HPP
#define SIMOL_INPUT_HPP

#include<yaml-cpp/yaml.h>
#include <getopt.h>
#include "Vector.hpp"

#include "io/CommandLine.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include "tools.hpp"

namespace simol{
  
  class Input{
    YAML::Node data;
		string inputPath_;
		ifstream inputFlux_;
    double positionMin_, positionMax_;
		ifstream inputSettings_;
		vector<dvec> initialPositions_, initialMomenta_;
  public:
    Input(CommandLine cmd);
		const string& inputPath() const;
		const std::ifstream& inputFlux() const;
		
		string simuTypeName() const;
		string parametersName() const;
    string outputFolderName() const;
		
		// Geometry
    int dimension() const;
    double length() const;
		
    // Mesh/Time
    double timeStepMin() const;
    double timeStepMax() const;
    double timeStep(size_t indexOfReplica=0) const;
		size_t nbOfIterations(size_t indexOfReplica=0) const;
		size_t nbOfThermalIterations(size_t indexOfReplica=0) const;
		size_t nbOfBurningIterations(size_t indexOfReplica=0) const;
		
		// Mesh/Replica
		size_t nbOfReplicas() const;
    
		// Physics/System
		string systemName() const;
    size_t nbOfParticles() const;
    double mass() const;
		bool doFileSettings() const;
		string settingsPath() const;
		dvec initialPosition(int const& i=0) const;
    dvec initialMomentum(int const& i=0) const;
		
		// Physics/Potential
    string potentialName() const;    
    //Sinus
    double amplitude() const;    
    //DoubleWell
    double height() const;
    double interWell() const;    
    //Harmonic
    double potentialStiffness() const;
    double potentialAlpha() const;
    double potentialBeta() const;
		
		// Physics/Model	
    string dynamicsName() const;
    double gamma() const;   
    double temperature(size_t indexOfReplica=0) const;
    double temperatureLeft(size_t indexOfReplica=0) const;
    double temperatureRight(size_t indexOfReplica=0) const;
    double beta(size_t indexOfReplica=0) const;
    double betaLeft(size_t indexOfReplica=0) const;
    double betaRight(size_t indexOfReplica=0) const;
    double externalForceMin() const;
    double externalForceMax() const;
    double externalForce(size_t indexOfReplica=0) const;
    double tauBending() const;
		double xi() const;
		int seed() const;
		double eta() const;
    
		
		
		// Output    
    size_t decorrelationNbOfIterations(size_t indexOfReplica=0) const;
    double decorrelationTime(size_t indexOfReplica=0) const;
		size_t outputPeriodNbOfIterations(size_t indexOfReplica=0) const;
    double outputPeriodTime(size_t indexOfReplica=0) const;
		size_t outputProfilePeriodNbOfIterations(size_t indexOfReplica=0) const;
    double outputProfilePeriodTime(size_t indexOfReplica=0) const;
		int nbOfAutocoPts() const;
		bool doFinalFlow() const;
		bool doFinalVelocity() const;
    
		// Controle Variate
    string controlVariateName() const;		
		string controlVariateCoeffsPath() const;
		
		//Galerkin
		bool doGalerkinCV() const;
		bool isGalerkin() const;
		string galerkinElts() const;
		size_t nbOfFourier() const;
		size_t nbOfHermite() const;
		
  };
  

}

#endif