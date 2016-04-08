#ifndef SIMOL_STATISTICS_HPP
#define SIMOL_STATISTICS_HPP

#include "Tools.hpp"
#include "Input.hpp"

namespace simol
{
  //template <typename T>
  class Statistics
  {
    DenseMatrix<double> sumValues_;
    DenseMatrix<double> lastValue_;
    DenseMatrix<int> nbValues_;
  public:
    Statistics(int nbRows = 1, int nbCols = 1);
    //virtual ~Statistics(){};
    void append(double value, int i = 0, int j = 0);
    double mean(int i = 0, int j = 0) const;
    Vector<double> meanVec(int i = 0) const;
    DenseMatrix<double> meanMat() const;
    const int& nbValues(int i = 0, int j = 0) const;
    const Vector<int> nbValuesVec(int i = 0) const;
    const DenseMatrix<int>& nbValuesMat() const;
    const double& lastValue(int i = 0, int j = 0) const;
    const Vector<double> lastValueVec(int i = 0) const;
    const DenseMatrix<double>& lastValueMat() const;
  };

  //
  //Calcule le profil de corrélation de deux observables A et B : <\Linv A, B>
  //Dans le cas où A = B on calcule la variance asymptotique de A
  //template <typename T>
  class AutocorrelationStats
  {
    int decorrelationNbOfIterations_;
    double decorrelationTime_;
		int nbOfAutocoPts_;
    int nbOfObservables_;
    Statistics statisticsValues_;        //calcule la moyenne de l'observable A
    Statistics statisticsRefValues_;     //calcule la moyenne de l'observable B
    Statistics statisticsMeanCorrelation_;   //calcule la moyenne de tous les A(x_i)B(x_0)
    Statistics statisticsCorrelation_;   //calcule la moyenne de A(t)B(0) pour t < decorrelationTime_
    double indexRef_;
  public:
		AutocorrelationStats();
    AutocorrelationStats(int decorrelationNbOfIterations, double decorrelationTime, int nbOfAutocoPts0, int nbOfObservables = 1);
    void append(double const& newValue, int iOfIteration, int iOfObservable = 0);
    void append(double const& newValue, int iOfIteration, int iOfObservable, double const& newRefValue);
    double operator()(int iOfIteration, int iOfObservable=0) const;
    const double& lastValue(int iOfObservable = 0) const;
    double mean(int iOfObservable = 0) const;
    int statisticsMeanCorrelationNbValues(int iOfObservable = 0) const;
    double integratedAutocorrelation(int iOfObservable = 0) const;
    Vector<double> integratedAutocorrelationVec() const;
    double integratedAutocorrelationUnbiased(int iOfObservable = 0) const;
    double variance(int iOfObservable = 0) const;
    double standardDeviation(int iOfObservable = 0) const;
  };
}




#endif
