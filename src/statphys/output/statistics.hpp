#ifndef SIMOL_STATISTICS_HPP
#define SIMOL_STATISTICS_HPP

#include "tools.hpp"
#include "input.hpp"

namespace simol
{
  //template <typename T>
  class Statistics
  {
    DenseMatrix<double> sumValues_;
    DenseMatrix<double> lastValue_;
    DenseMatrix<size_t> nbValues_;
  public:
    Statistics(size_t nbRows = 1, size_t nbCols = 1);
    //virtual ~Statistics(){};
    void append(double value, size_t i = 0, size_t j = 0);
    double mean(size_t i = 0, size_t j = 0) const;
    Vector<double> meanVec(size_t i = 0) const;
    DenseMatrix<double> meanMat() const;
    const size_t& nbValues(size_t i = 0, size_t j = 0) const;
    const Vector<size_t> nbValuesVec(size_t i = 0) const;
    const DenseMatrix<size_t>& nbValuesMat() const;
    const double& lastValue(size_t i = 0, size_t j = 0) const;
    const Vector<double> lastValueVec(size_t i = 0) const;
    const DenseMatrix<double>& lastValueMat() const;
  };

  //
  //Calcule le profil de corrélation de deux observables A et B : <\Linv A, B>
  //Dans le cas où A = B on calcule la variance asymptotique de A
  //template <typename T>
  class AutocorrelationStats
  {
    size_t decorrelationNbOfIterations_;
    double decorrelationTime_;
		int nbOfAutocoPts_;
    size_t nbOfObservables_;
    Statistics statisticsValues_;        //calcule la moyenne de l'observable A
    Statistics statisticsRefValues_;     //calcule la moyenne de l'observable B
    Statistics statisticsMeanCorrelation_;   //calcule la moyenne de tous les A(x_i)B(x_0)
    Statistics statisticsCorrelation_;   //calcule la moyenne de A(t)B(0) pour t < decorrelationTime_
    double indexRef_;
  public:
		AutocorrelationStats();
    AutocorrelationStats(size_t decorrelationNbOfIterations, double decorrelationTime, int nbOfAutocoPts0, size_t nbOfObservables = 1);
    void append(double const& newValue, size_t iOfIteration, size_t iOfObservable = 0);
    void append(double const& newValue, size_t iOfIteration, size_t iOfObservable, double const& newRefValue);
    double operator()(size_t iOfIteration, size_t iOfObservable=0) const;
    const double& lastValue(size_t iOfObservable = 0) const;
    double mean(size_t iOfObservable = 0) const;
    size_t statisticsMeanCorrelationNbValues(size_t iOfObservable = 0) const;
    double integratedAutocorrelation(size_t iOfObservable = 0) const;
    Vector<double> integratedAutocorrelationVec() const;
    double integratedAutocorrelationUnbiased(size_t iOfObservable = 0) const;
    double variance(size_t iOfObservable = 0) const;
    double standardDeviation(size_t iOfObservable = 0) const;
  };
}




#endif
