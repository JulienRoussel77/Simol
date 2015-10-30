#pragma once

#include "input.hpp"
#include <iostream>
#include <vector>
using std::vector;
using std::cout; 
using std::endl;


namespace simol
{
  template <typename T>
  class Statistics
  {
    vector<T> sumValues_;
    vector<T> lastValue_;
    vector<int> nbValues_;
  public:
    Statistics(int nbIndices = 1);
    //virtual ~Statistics(){};
    void append(T value, int i = 0);
    T mean(int i = 0) const;
    const int& nbValues(int i = 0) const;
    const T& lastValue(int i = 0) const;
  };
  
  template <typename T>
  class AutocorrelationStats
  {
    size_t decorrelationNumberOfIterations_;
    //double timeStep_;
    Statistics<T> statisticsValues_;
    Statistics<T> statisticsRefValues_;
    Statistics<double> statisticsMeanCorrelation_; 
    Statistics<double> statisticsMeanCorrelationTemp_; 
    Statistics<double> statisticsCorrelation_; 
    double indexRef_;
    //T valueRef_;
    //T lastValue_;
  public:
    AutocorrelationStats(size_t decorrelationNumberOfIterations);
    //const double& timeStep() const;
    void append(T const& newValue, int indexOfIteration);
    double operator()(int i) const;
    const T& lastValue(int i = 0) const;
    T mean() const;
    size_t statisticsMeanCorrelationNbValues(int i = 0) const;
    double integratedAutocorrelation() const;
    double integratedAutocorrelationUnbiased() const;
    double variance() const;
    double standardDeviation() const;
  };
}


namespace simol
{
  template <class T>
  Statistics<T>::Statistics(int nbIndices):sumValues_(nbIndices, 0), lastValue_(nbIndices, 0), nbValues_(nbIndices, 0){};
  
  template <class T>
  void Statistics<T>::append(T value, int i)
  {
    lastValue_[i] = value;
    if (nbValues_[i] == 0)
      sumValues_[i] = value;
    else
      sumValues_[i] += value;
    nbValues_[i]++;
  }

  template <class T>  
  T Statistics<T>::mean(int i) const
  {
    if (nbValues_[i] == 0)
      return 0;
    else
      return sumValues_[i] / nbValues_[i];
  }
  
  template <class T>
  const int& Statistics<T>::nbValues(int i) const
  {
    return nbValues_[i];
  }
  
  template <class T>
  const T& Statistics<T>::lastValue(int i) const
  {
    return lastValue_[i];
  }
  
  
  
  
  template <class T>
  AutocorrelationStats<T>::AutocorrelationStats(size_t decorrelationNumberOfIterations):
    decorrelationNumberOfIterations_(decorrelationNumberOfIterations),
    //timeStep_(timeStep),
    statisticsValues_(),
    statisticsRefValues_(),
    statisticsMeanCorrelation_(),
    statisticsMeanCorrelationTemp_(),
    statisticsCorrelation_(decorrelationNumberOfIterations_),
    indexRef_(0)
  {}
  
  template <class T>
  void AutocorrelationStats<T>::append(T const& newValue, int indexOfIteration)
  {
    cout << "AutocorrelationStats<T>::append(T newValue, int indexOfIteration) not implemented !" << endl;
  
    exit(1);
  }
  
  template <class T>
  double AutocorrelationStats<T>::operator()(int i) const
  {
    return statisticsCorrelation_.mean(i);
  }
  
  template <class T>
  const T& AutocorrelationStats<T>::lastValue(int i) const
  {
    return statisticsValues_.lastValue(i);
  }
  
  template <class T>
  T AutocorrelationStats<T>::mean() const
  {
    return statisticsValues_.mean();
  }
  
  template <class T>
  size_t AutocorrelationStats<T>::statisticsMeanCorrelationNbValues(int i) const
  {
    return statisticsMeanCorrelation_.nbValues(i);
  }
  
  template <class T>
  double AutocorrelationStats<T>::integratedAutocorrelation() const
  {
    return statisticsMeanCorrelation_.mean();// * decorrelationNumberOfIterations_ * timeStep();
  }
  
  template <class T>
  double AutocorrelationStats<T>::integratedAutocorrelationUnbiased() const
  {
    cout << "AutocorrelationStats<T>::integratedAutocorrelationUnbiased() not implemented !" << endl;
    exit(1);
    return 0;
  }
  
  template <class T>
  double AutocorrelationStats<T>::variance() const
  {
    return 2 * integratedAutocorrelationUnbiased();
  }
  
    template <class T>
  double AutocorrelationStats<T>::standardDeviation() const
  {
    return sqrt(variance());
  }

  
}