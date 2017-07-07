#ifndef SIMOL_STATISTICS_HPP
#define SIMOL_STATISTICS_HPP

#include "simol/statphys/Tools.hpp"
#include "simol/statphys/input/Input.hpp"

namespace simol
{
  ///Compute on the line the mean of a set of iid variables
  ///The class is overloaded so that it can deal with nbRows x nbCols observables simultaneously
  class Statistics
  {
  protected:
    DMat sumValues_;
    DMat lastValue_;
    DMatInt nbValues_;
    DMat iidVar_;
  public:
    Statistics(int nbRows = 1, int nbCols = 1);
    //virtual ~Statistics(){};
    int nbOfRows() const;
    int nbOfCols() const;
    virtual void append(double value, int i = 0, int j = 0);
    double mean(int i = 0, int j = 0) const;
    DVec meanVec(int i = 0) const;
    DMat meanMat() const;
    const long int& nbValues(int i = 0, int j = 0) const;
    const DVecInt nbValuesVec(int i = 0) const;
    const DMatInt& nbValuesMat() const;
    const double& lastValue(int i = 0, int j = 0) const;
    const DVec lastValueVec(int i = 0) const;
    const DMat& lastValueMat() const;      
  };
  
  ///Compute on the line the mean and variance of a set of iid variables
  ///The variance is computed using a recursive formula so only a few scalar need to be stocked  
  class IIDStats : public Statistics
  {
    DMat sumVar_; // We store a cumulated version of the variance to avoid numerical instabilities
  public:
    IIDStats(int nbRows = 1, int nbCols = 1);
    void append(double value, int i = 0, int j = 0);
    double variance(int i = 0, int j = 0 ) const;
    double stdDev(int i = 0, int j = 0) const;
  };
  
  //
  //Calcule le profil de corrélation de deux observables A et B : <\Linv A, B>
  //Dans le cas où A = B on calcule la variance asymptotique de A
  //template <typename T>
  class CorrelationStats
  {
  protected: 
    int decorrelationNbOfSteps_;
    double timeStep_;
    int nbOfAutocoPts_;
    int nbOfObservables_;
    Statistics statsValues_;            //calcule la moyenne de l'observable A
    IIDStats statsRefValues_;             //calcule la moyenne de l'observable B
    IIDStats statsIntegratedCorrelation_;   //calcule la moyenne de tous les A(x_i)B(x_0)
    IIDStats statsCorrelation_;      //calcule la moyenne de A(t)B(0) pour t < decorrelationTime_
    long int indexRef_;                      //index of reference for B(0)
    double currentCorrelationIntegral_;      //contains the correlations cumulated since the last index of reference
  public:
    CorrelationStats();
    CorrelationStats(int decorrelationNbOfSteps, double timeStep0, int nbOfAutocoPts0, int nbOfObservables = 1);
    const int& decorrelationNbOfSteps() const;
    const double& timeStep() const;
    double decorrelationTime() const;
    //void append(double const& newValue, long int iOfStep, int iOfObservable = 0);
    virtual void append(double const& newValue, long int iOfStep, int iOfObservable, double const& newRefValue);
    //double operator()(long int iOfStep, int iOfObservable = 0) const;
    virtual double correlationAtSpan(long int iOfSpan, int iOfObservable = 0) const;
    virtual double unbiasedCorrelationAtSpan(long int iOfSpan, int iOfObservable = 0) const;
    const double& lastValue(int iOfObservable = 0) const;
    double mean(int iOfObservable = 0) const;
    const long int& nbValues(int i = 0, int j = 0) const;
    int statsIntegratedCorrelationNbValues(int iOfObservable = 0) const;
    double integratedCorrelation(int iOfObservable = 0) const;
    DVec integratedCorrelationVec() const;
    virtual double integratedCorrelationUnbiased(int iOfObservable = 0) const;
    double varIntegratedCorrelation(int iOfObservable = 0) const;
    double varCorrelationAtSpan(int iOfSpan, int iOfObservable = 0) const;
    double stdDevCorrelationAtSpan(int iOfSpan, int iOfObservable = 0) const;
    double stdErrorCorrelationAtSpan(int iOfSpan, int iOfObservable = 0) const;
  };
  
  ///
  ///Case when A = B: typically an asymptotic variance estimation
  class AutocorrelationStats : public CorrelationStats
  {
  public:
    AutocorrelationStats();
    AutocorrelationStats(int decorrelationNbOfSteps, double timeStep0, int nbOfAutocoPts0, int nbOfObservables = 1);
    void append(double const& newValue, long int iOfStep, int iOfObservable=0);
    double unbiasedCorrelationAtSpan(long int iOfStep, int iOfObservable = 0) const;
    double integratedCorrelationUnbiased(int iOfObservable=0) const;
    double variance(int iOfObservable = 0) const;
    double stdDev(int iOfObservable = 0) const;
    double varOfVar(int iOfObservable = 0) const;
    double stdDevOfVar(int iOfObservable = 0) const;
    //double stdErrorOfVar(int iOfObservable = 0) const;
  };
}




#endif
