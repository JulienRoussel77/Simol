#ifndef SIMOL_OBSERVABLE_HPP
#define SIMOL_OBSERVABLE_HPP

#include <iomanip>
using std::setw;
#include "simol/statphys/Tools.hpp"
#include "Statistics.hpp"
#include "simol/statphys/system/Particle.hpp"

namespace simol
{
  //class Observable;
  //Observable* createObservable(const Input& input);
  
 class Observable
  {
  //friend Observable* createObservable(const Input& input);
  public:
      int decorrelationNbOfSteps_;
      double timeStep_;
      int printPeriodNbOfSteps_;
      int nbOfAutocoPts_;

      double currentValue_;
      AutocorrelationStats autocoStats_;
      
      ofstream outFlux_;      
      ofstream* outFluxCorrelation_;
    public:
      Observable(Input const& input, const string& outPath);
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
      
      /*void computeGeneratorHamiltonian(vector<Particle> const& configuration){};
      void computeGeneratorOverdamped(vector<Particle> const& configuration, double beta){};
      void computeGeneratorLangevin(vector<Particle> const& configuration, double beta, double gamma){};
      void computeGeneratorBoundarylangevin(vector<Particle> const& configuration, double betaLeft, double betaRight, double gamma){};*/
  
      virtual void updateHamiltonian(vector<Particle> const& /*configuration*/){};
      virtual void updateOverdamped(vector<Particle> const& /*configuration*/, double /*beta*/){};
      virtual void updateLangevin(vector<Particle> const& /*configuration*/, double /*beta*/, double /*gamma*/){};
      virtual void updateBoundaryLangevin(vector<Particle> const& /*configuration*/, double /*betaLeft*/, double /*betaRight*/, double /*gamma*/){};

      bool doOutput(long int iOfStep) const;
      //bool doFinal() const;
      virtual void display(long int iOfStep);
      //virtual void finalDisplay(vector<double> parameters);
      void displayCorrelations(long int iOfStep);
  };
}

#endif