#pragma once

#include "statistics.hpp"
#include "potential.hpp"
#include <iostream>
#include <fstream>
 
namespace simol
{
  
  class ControlVariate;
  ControlVariate* createControlVariate(Input const& input, Potential& potential, size_t indexOfReplica=1);
  
  class ControlVariate
  {
  friend ControlVariate* createControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
  protected:
    AutocorrelationStats<double> statsB_;
    AutocorrelationStats<double> statsD_;
    AutocorrelationStats<double> statsObservable_;
    AutocorrelationStats<double> statsBetterObservable_;
    
    AutocorrelationStats<double> statsGeneratorOnBasis_;
    /*double lastValueB1_;
    double estimatorB1_;
    int nbValuesB1_;*/
    //Statistics<double> estimatorB2_;

    /*double lastValueD_;    
    double estimatorD_;
    int nbValuesD_;
    double lastValueObservable_;    
    double meanObservable_;
    int nbValuesObservable_;
    double lastValueBetterObservable_;    
    double meanBetterObservable_;
    int nbValuesBetterObservable_;*/
    
    size_t decorrelationNumberOfIterations_;
    Potential* potential_;
  public:
    ControlVariate(Input const& input, Potential& potential);
    double potential(dvec const& position) const;
    dvec potentialDerivative(dvec const& position) const;
    double potentialLaplacian(dvec const& position) const;
    virtual double basisFunction(dvec const& position, dvec const& momentum) const = 0;
    //virtual double generatorOnBasisFunction(dvec const& position, dvec const& momentum) const = 0;
    
    //virtual double observableFunction(dvec const& position, dvec const& momentum) = 0;
    void appendToB(double observable, dvec const& position, dvec const& momentum, size_t indexOfIteration);
    void appendToD(double generatorOnBasisFunction, dvec const& position, dvec const& momentum, size_t indexOfIteration);
    void appendToObservable(double observable, size_t indexOfIteration);
    void appendToBetterObservable(double observable, double generatorOnBasisFunction, dvec const& position, dvec const& momentum, size_t indexOfIteration);
    virtual void update(double observable, double generatorOnBasisFunction, dvec const& position, dvec const& momentum, size_t indexOfIteration);
    
    virtual double lastB() const;
    virtual double meanB() const;
    virtual double stdDeviationB() const;
    virtual double lastD() const;
    virtual double meanD() const;
    virtual double stdDeviationD() const;
    virtual double lastObservable() const;
    virtual double meanObservable() const;
    virtual double stdDeviationObservable() const;
    virtual double lastBetterObservable() const;
    virtual double meanBetterObservable() const;
    virtual double stdDeviationBetterObservable() const;
    virtual double lastGeneratorOnBasis() const;
    virtual double meanGeneratorOnBasis() const;
    virtual double stdDeviationGeneratorOnBasis() const;
    
    virtual double autocorrelation(size_t indexDifference) const;
    
    virtual double laplacienQ(dvec const& position, dvec const& momentum) const = 0;
    virtual dvec gradientQ(dvec const& position, dvec const& momentum) const = 0;
    virtual double laplacienP(dvec const& position, dvec const& momentum) const = 0;
    virtual dvec gradientP(dvec const& position, dvec const& momentum) const = 0;
  };
  
  
  class SinusControlVariate : public ControlVariate
  {
  public:
    SinusControlVariate(Input const& input, Potential& potential);
    double basisFunction(dvec const& position, dvec const& momentum) const;
    //double generatorOnBasisFunction(dvec const& position, dvec const& momentum) const;
    virtual double laplacienQ(dvec const& position, dvec const& momentum) const;
    virtual dvec gradientQ(dvec const& position, dvec const& momentum) const;
    virtual double laplacienP(dvec const& position, dvec const& momentum) const;
    virtual dvec gradientP(dvec const& position, dvec const& momentum) const;
  };
  
  class SinExpControlVariate : public ControlVariate
  {
  public:
    SinExpControlVariate(Input const& input, Potential& potential);
    double basisFunction(dvec const& position, dvec const& momentum) const;
    //double generatorOnBasisFunction(dvec const& position, dvec const& momentum) const;
    virtual double laplacienQ(dvec const& position, dvec const& momentum) const;
    virtual dvec gradientQ(dvec const& position, dvec const& momentum) const;
    virtual double laplacienP(dvec const& position, dvec const& momentum) const;
    virtual dvec gradientP(dvec const& position, dvec const& momentum) const;
  };
  
  class NoControlVariate : public ControlVariate
  {
  public:
    NoControlVariate(Input const& input, Potential& potential);
    double basisFunction(dvec const& position, dvec const& momentum) const;
    //double generatorOnBasisFunction(dvec const& position, dvec const& momentum) const;
    virtual double laplacienQ(dvec const& position, dvec const& momentum) const;
    virtual dvec gradientQ(dvec const& position, dvec const& momentum) const;
    virtual double laplacienP(dvec const& position, dvec const& momentum) const;
    virtual dvec gradientP(dvec const& position, dvec const& momentum) const;
    void update(double observable, double generatorOnBasisFunction, dvec const& position, dvec const& momentum, size_t indexOfIteration);

  };
  
  
}