#ifndef SIMOL_OBSERVABLE_HPP
#define SIMOL_OBSERVABLE_HPP

#include <iomanip>
using std::setw;
#include "simol/statphys/Tools.hpp"
#include "Statistics.hpp"
#include "simol/statphys/system/Particle.hpp"
#include "simol/statphys/system/System.hpp"

namespace simol
{
  
 class Observable
  {
  //friend Observable* createObservable(const Input& input);
  public:
      int decorrelationNbOfSteps_;
      double timeStep_;
      int printPeriodNbOfSteps_;
      int nbOfAutocoPts_;

      //Contains the value that is being computed, set to zero at each iteration
      double currentValue_;
      AutocorrelationStats autocoStats_;
      
      string outPath_;
      ofstream outFlux_;      
      shared_ptr<ofstream> outFluxCorrelation_;
    public:
      Observable(Input const& input, int idObs);
      virtual ~Observable();
      
      const double& timeStep() const;
      const int& decorrelationNbOfSteps() const;
      double decorrelationTime() const;      
      int const& nbOfAutocoPts() const;
      double autocoPtsPeriod() const;
      double printPeriodTime() const;
      const int& printPeriodNbOfSteps() const;
      double time() const;
      
      virtual int nbOfFourier() const
      {return 0;}
      
      virtual int nbOfHermite() const
      {return 0;}
      
      string const& outPath() const;
      ofstream& outFlux();
      ofstream& outFluxCorrelation();
      
      bool doComputeCorrelations() const;

      virtual const double& currentValue() const;
      virtual double& currentValue();
      virtual double lastValue() const;
      virtual double mean() const;
      virtual double variance() const;
      virtual double stdDev() const;
      virtual double varOfVar() const;
      virtual double stdDevOfVar() const;

      virtual double correlationAtSpan(int iOfSpan) const;
      virtual double unbiasedCorrelationAtSpan(int iOfSpan) const;
      virtual double varCorrelationAtSpan(int iOfSpan) const;
      virtual double stdDevCorrelationAtSpan(int iOfSpan) const;

      virtual void append(double value, long int iOfStep);
      virtual void appendCurrent(long int iOfStep);
      
      /*void computeGeneratorHamiltonian(){};
      void computeGeneratorOverdamped(System const& syst, double beta){};
      void computeGeneratorLangevin(System const& syst, double beta, double gamma){};
      void computeGeneratorBoundarylangevin(System const& syst, double betaLeft, double betaRight, double gamma){};*/
  
      virtual void updateHamiltonian(System const&){};
      virtual void updateOverdamped(System const&, double /*beta*/){};
      virtual void updateLangevin(System const&, double /*beta*/, double /*gamma*/){};
      virtual void updateBoundaryLangevin(System const&, double /*betaLeft*/, double /*betaRight*/, double /*gamma*/){};

      virtual bool doOutput(long int iOfStep) const;
      //bool doFinal() const;
      virtual void display(long int iOfStep);
      virtual void displayFinalValues(ofstream& out);

      virtual void displayCorrelations(long int iOfStep);
  };
}

#endif