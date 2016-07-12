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
    DenseMatrix<double> sumValues_;
    DenseMatrix<double> lastValue_;
    DenseMatrix<long int> nbValues_;
    DenseMatrix<double> iidVar_;
  public:
    Statistics(int nbRows = 1, int nbCols = 1);
    //virtual ~Statistics(){};
    virtual void append(double value, int i = 0, int j = 0);
    double mean(int i = 0, int j = 0) const;
    Vector<double> meanVec(int i = 0) const;
    DenseMatrix<double> meanMat() const;
    const long int& nbValues(int i = 0, int j = 0) const;
    const Vector<long int> nbValuesVec(int i = 0) const;
    const DenseMatrix<long int>& nbValuesMat() const;
    const double& lastValue(int i = 0, int j = 0) const;
    const Vector<double> lastValueVec(int i = 0) const;
    const DenseMatrix<double>& lastValueMat() const;      
  };
  
  ///Compute on the line the mean and variance of a set of iid variables
  ///The variance is computed using a recursive formula so only a few scalar need to be stocked  
  class IIDStats : public Statistics
  {
    DenseMatrix<double> iidVar_;
  public:
    IIDStats(int nbRows = 1, int nbCols = 1);
    void append(double value, int i = 0, int j = 0);
    const double& variance(int i = 0, int j = 0 ) const;
    double stdDeviation(int i = 0, int j = 0) const;
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
    void append(double const& newValue, long int iOfStep, int iOfObservable, double const& newRefValue);
    double operator()(long int iOfStep, int iOfObservable = 0) const;
    const double& lastValue(int iOfObservable = 0) const;
    double mean(int iOfObservable = 0) const;
    int statsIntegratedCorrelationNbValues(int iOfObservable = 0) const;
    double integratedCorrelation(int iOfObservable = 0) const;
    Vector<double> integratedCorrelationVec() const;
    double integratedCorrelationUnbiased(int iOfObservable = 0) const;
    double varIntegratedCorrelation(int iOfObservable = 0) const;
    double varCorrelation(int indexDifference, int iOfObservable = 0) const;
    double stdDeviationCorrelation(int indexDifference, int iOfObservable = 0) const;
    double stdErrorCorrelation(int indexDifference, int iOfObservable = 0) const;
  };
  
  ///
  ///Case when A = B: typically an asymptotic variance estimation
  class AutocorrelationStats : public CorrelationStats
  {
  public:
    AutocorrelationStats();
    AutocorrelationStats(int decorrelationNbOfSteps, double timeStep0, int nbOfAutocoPts0, int nbOfObservables = 1);
    void append(double const& newValue, long int iOfStep, int iOfObservable=0);
    double variance(int iOfObservable = 0) const;
    double standardDeviation(int iOfObservable = 0) const;
    double varianceOfVariance(int iOfObservable = 0) const;
    double stdDevOfVariance(int iOfObservable = 0) const;
    double stdErrorOfVariance(int iOfObservable = 0) const;
  };
}




#endif
