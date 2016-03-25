#ifndef SIMOL_STATISTICS_HPP
#define SIMOL_STATISTICS_HPP

#include "tools.hpp"
#include "input.hpp"

namespace simol
{
  //template <typename T>
  class Statistics
  {
    /*vector sumValues_;
    vector lastValue_;
    vector<int> nbValues_;*/

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
    const Vector<size_t>& nbValuesVec(size_t i = 0) const;
    const DenseMatrix<size_t>& nbValuesMat() const;
    const double& lastValue(size_t i = 0, size_t j = 0) const;
    const Vector<double>& lastValueVec(size_t i = 0) const;
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


namespace simol
{
  Statistics::Statistics(size_t nbRows, size_t nbCols):sumValues_(nbRows, nbCols), lastValue_(nbRows, nbCols), nbValues_(nbRows, nbCols)
  {
		//cout << "Statistics " << nbRows << " X " << nbCols << endl;
    sumValues_.fill(0);
    lastValue_.fill(0);
    nbValues_.fill(0);
  };

  void Statistics::append(double value, size_t i, size_t j)
  {
    lastValue_(i,j) = value;
    if (nbValues_(i,j) == 0)
      sumValues_(i,j) = value;
    else
      sumValues_(i,j) += value;
    nbValues_(i,j)++;
  }


  double Statistics::mean(size_t i, size_t j) const
  {
    if (nbValues_(i,j) == 0)
      return 0;
    else 
      return sumValues_(i,j) / nbValues(i,j);
  }
  
  Vector<double> Statistics::meanVec(size_t i) const
  {
    //return (sumValues_.array() / nbValues_.array().cast<double>()).matrix();
    return piecewiseDivision(sumValues_.column(i), nbValues_.column(i));
  }
  
  DenseMatrix<double> Statistics::meanMat() const
  {
    //return (sumValues_.array() / nbValues_.array().cast<double>()).matrix();
    return piecewiseDivision(sumValues_, nbValues_);
  }


  const size_t& Statistics::nbValues(size_t i, size_t j) const
  {
    return nbValues_(i,j);
  }

  const Vector<size_t>& Statistics::nbValuesVec(size_t i) const
  {
    return nbValues_.column(i);
  }

  const DenseMatrix<size_t>& Statistics::nbValuesMat() const
  {
    return nbValues_;
  }


  const double& Statistics::lastValue(size_t i, size_t j) const
  {
    return lastValue_(i,j);
  }

  const Vector<double>& Statistics::lastValueVec(size_t i) const
  {
    return lastValue_.column(i);
  }

  const DenseMatrix<double>& Statistics::lastValueMat() const
  {
    return lastValue_;
  }




  AutocorrelationStats::AutocorrelationStats():
    decorrelationNbOfIterations_(0),
    decorrelationTime_(0),
    nbOfAutocoPts_(0),
    nbOfObservables_(0),
    //timeStep_(timeStep),
    statisticsValues_(0),
    statisticsRefValues_(0),
    statisticsMeanCorrelation_(0),
    //statisticsMeanCorrelationTemp_(),
    statisticsCorrelation_(0, 0),
    indexRef_(0)
  {//cout << "AutocorrelationStats()" << endl;
	}


  AutocorrelationStats::AutocorrelationStats(size_t decorrelationNbOfIterations0, double decorrelationTime, int nbOfAutocoPts0, size_t nbOfObservables):
    decorrelationNbOfIterations_(decorrelationNbOfIterations0),
    decorrelationTime_(decorrelationTime),
    nbOfAutocoPts_(nbOfAutocoPts0),
    nbOfObservables_(nbOfObservables),
    //timeStep_(timeStep),
    statisticsValues_(nbOfObservables),
    statisticsRefValues_(nbOfObservables),
    statisticsMeanCorrelation_(nbOfObservables),
    //statisticsMeanCorrelationTemp_(),
    statisticsCorrelation_(nbOfAutocoPts_, nbOfObservables),
    indexRef_(0)
  {//cout << "AutocorrelationStats(size_t decorrelation : nbObs = " << decorrelationTime << endl;
	}


  void AutocorrelationStats::append(double const& newValue, size_t iOfIteration, size_t iOfObservable)
  {
    append(newValue, iOfIteration, iOfObservable, newValue);
  }


  double AutocorrelationStats::operator()(size_t iOfIteration, size_t iOfObservable) const
  {
    return statisticsCorrelation_.mean(iOfIteration, iOfObservable);
  }


  const double& AutocorrelationStats::lastValue(size_t iOfObservable) const
  {
    return statisticsValues_.lastValue(iOfObservable);
  }


  double AutocorrelationStats::mean(size_t iOfObservable) const
  {
    return statisticsValues_.mean(iOfObservable);
  }


  size_t AutocorrelationStats::statisticsMeanCorrelationNbValues(size_t iOfObservable) const
  {
    return statisticsMeanCorrelation_.nbValues(iOfObservable);
  }


  double AutocorrelationStats::integratedAutocorrelation(size_t iOfObservable) const
  {
    return statisticsMeanCorrelation_.mean(iOfObservable) * decorrelationTime_;
  }


  Vector<double> AutocorrelationStats::integratedAutocorrelationVec() const
  {
    return statisticsMeanCorrelation_.meanMat().column(0) * decorrelationTime_;
  }




  double AutocorrelationStats::variance(size_t iOfObservable) const
  {
    return 2 * integratedAutocorrelationUnbiased(iOfObservable);
  }

  
  double AutocorrelationStats::standardDeviation(size_t iOfObservable) const
  {
    return sqrt(variance(iOfObservable));
  }


}

#endif
