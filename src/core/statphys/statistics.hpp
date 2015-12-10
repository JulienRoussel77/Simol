#pragma once

#include "input.hpp"
#include <iostream>
#include <vector>
using std::vector;
using std::cout; 
using std::endl;
using Eigen::Matrix;
using Eigen::Dynamic;

namespace simol
{
  template <typename T>
  class Statistics
  {
    /*vector<T> sumValues_;
    vector<T> lastValue_;
    vector<int> nbValues_;*/
    
    Matrix<T, Dynamic, Dynamic> sumValues_;
    Matrix<T, Dynamic, Dynamic> lastValue_;
    Matrix<size_t, Dynamic, Dynamic> nbValues_;    
  public:
    Statistics(size_t nbRows = 1, size_t nbCols = 1);
    //virtual ~Statistics(){};
    void append(T value, size_t i = 0, size_t j = 0);
    T mean(size_t i = 0, size_t j = 0) const;
    Matrix<T, Dynamic, Dynamic> meanMat() const;
    const size_t& nbValues(size_t i = 0, size_t j = 0) const;
    const Matrix<size_t, Dynamic, Dynamic>& nbValuesMat() const;
    const T& lastValue(size_t i = 0, size_t j = 0) const;
    const Matrix<T, Dynamic, Dynamic>& lastValueMat() const;
  };
  
  //
  //Calcule le profil de corrélation de deux observables A et B : <\Linv A, B>
  //Dans le cas où A = B on calcule la variance asymptotique de A
  template <typename T>
  class AutocorrelationStats
  {
    size_t decorrelationNumberOfIterations_;
    double decorrelationTime_;
    size_t nbOfObservables_;
    //double timeStep_;
    Statistics<T> statisticsValues_;        //calcule la moyenne de l'observable A
    Statistics<T> statisticsRefValues_;     //calcule la moyenne de l'observable B
    Statistics<double> statisticsMeanCorrelation_;   //calcule la moyenne de tous les A(x_i)B(x_0)
    //Statistics<double> statisticsMeanCorrelationTemp_;  //calcule 1/n sum_i^n A(x_i)B(x_0) pour une trajectoire
    Statistics<double> statisticsCorrelation_;   //calcule la moyenne de A(t)B(0) pour t < decorrelationTime_
    double indexRef_;
    //T valueRef_;
    //T lastValue_;
  public:
    AutocorrelationStats(size_t decorrelationNumberOfIterations, double decorrelationTime, size_t nbOfObservables = 1);
    //const double& timeStep() const;
    void append(T const& newValue, size_t iOfIteration, size_t iOfObservable = 0);
    void append(T const& newValue, size_t iOfIteration, size_t iOfObservable, T const& newRefValue);
    double operator()(size_t iOfIteration, size_t iOfObservable=0) const;
    const T& lastValue(size_t iOfObservable = 0) const;
    T mean(size_t iOfObservable = 0) const;
    size_t statisticsMeanCorrelationNbValues(size_t iOfObservable = 0) const;
    double integratedAutocorrelation(size_t iOfObservable = 0) const;
    Eigen::MatrixXd integratedAutocorrelationMat() const;
    double integratedAutocorrelationUnbiased(size_t iOfObservable = 0) const;
    double variance(size_t iOfObservable = 0) const;
    double standardDeviation(size_t iOfObservable = 0) const;
  };
}


namespace simol
{
  template <class T>
  Statistics<T>::Statistics(size_t nbRows, size_t nbCols):sumValues_(nbRows, nbCols), lastValue_(nbRows, nbCols), nbValues_(nbRows, nbCols)
  {
    sumValues_.fill(0);
    lastValue_.fill(0);
    nbValues_.fill(0);
  };
  
  template <class T>
  void Statistics<T>::append(T value, size_t i, size_t j)
  {
    lastValue_(i,j) = value;
    if (nbValues_(i,j) == 0)
      sumValues_(i,j) = value;
    else
      sumValues_(i,j) += value;
    nbValues_(i,j)++;
  }

  template <class T>  
  T Statistics<T>::mean(size_t i, size_t j) const
  {
    if (nbValues_(i,j) == 0)
      return 0;
    else
      return sumValues_(i,j) / nbValues_(i,j);
  }
  
  template <class T>  
  Matrix<T, Dynamic, Dynamic> Statistics<T>::meanMat() const
  {
    return (sumValues_.array() / nbValues_.array().cast<double>()).matrix();
  }
  
  template <class T>
  const size_t& Statistics<T>::nbValues(size_t i, size_t j) const
  {
    return nbValues_(i,j);
  }
  
  template <class T>
  const Matrix<size_t, Dynamic, Dynamic>& Statistics<T>::nbValuesMat() const
  {
    return nbValues_;
  }
  
  template <class T>
  const T& Statistics<T>::lastValue(size_t i, size_t j) const
  {
    return lastValue_(i,j);
  }
  
  template <class T>
  const Matrix<T, Dynamic, Dynamic>& Statistics<T>::lastValueMat() const
  {
    return lastValue_;
  }
  
  
  
  
  template <class T>
  AutocorrelationStats<T>::AutocorrelationStats(size_t decorrelationNumberOfIterations, double decorrelationTime, size_t nbOfObservables):
    decorrelationNumberOfIterations_(decorrelationNumberOfIterations),
    decorrelationTime_(decorrelationTime),
    nbOfObservables_(nbOfObservables),
    //timeStep_(timeStep),
    statisticsValues_(nbOfObservables),
    statisticsRefValues_(nbOfObservables),
    statisticsMeanCorrelation_(nbOfObservables),
    //statisticsMeanCorrelationTemp_(),
    statisticsCorrelation_(decorrelationNumberOfIterations_, nbOfObservables),
    indexRef_(0)
  {}
  
  template <class T>
  void AutocorrelationStats<T>::append(T const& newValue, size_t iOfIteration, size_t iOfObservable)
  {
    append(newValue, iOfIteration, iOfObservable, newValue);
  }
  
  template <class T>
  void AutocorrelationStats<T>::append(T const& newValue, size_t iOfIteration, size_t iOfObservable, T const& newRefValue)
  {
    cout << "AutocorrelationStats<T>::append(T newValue, size_t iOfIteration) not implemented !" << endl;
  
    assert(false);
  }
  
  template <class T>
  double AutocorrelationStats<T>::operator()(size_t iOfIteration, size_t iOfObservable) const
  {
    return statisticsCorrelation_.mean(iOfIteration, iOfObservable);
  }
  
  template <class T>
  const T& AutocorrelationStats<T>::lastValue(size_t iOfObservable) const
  {
    return statisticsValues_.lastValue(iOfObservable);
  }
  
  template <class T>
  T AutocorrelationStats<T>::mean(size_t iOfObservable) const
  {
    return statisticsValues_.mean(iOfObservable);
  }
  
  template <class T>
  size_t AutocorrelationStats<T>::statisticsMeanCorrelationNbValues(size_t iOfObservable) const
  {
    return statisticsMeanCorrelation_.nbValues(iOfObservable);
  }
  
  template <class T>
  double AutocorrelationStats<T>::integratedAutocorrelation(size_t iOfObservable) const
  {
    return statisticsMeanCorrelation_.mean(iOfObservable) * decorrelationTime_;
  }
  
  template <class T>
  Eigen::MatrixXd AutocorrelationStats<T>::integratedAutocorrelationMat() const
  {
    return statisticsMeanCorrelation_.meanMat() * decorrelationTime_;
  }
  
  template <class T>
  double AutocorrelationStats<T>::integratedAutocorrelationUnbiased(size_t iOfObservable) const
  {
    cout << "AutocorrelationStats<T>::integratedAutocorrelationUnbiased() not implemented !" << endl;
    exit(1);
    return 0;
  }
  
  template <class T>
  double AutocorrelationStats<T>::variance(size_t iOfObservable) const
  {
    return 2 * integratedAutocorrelationUnbiased(iOfObservable);
  }
  
    template <class T>
  double AutocorrelationStats<T>::standardDeviation(size_t iOfObservable) const
  {
    return sqrt(variance(iOfObservable));
  }

  
}