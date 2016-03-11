#include "statistics.hpp"

using std::min;
using std::max;


namespace simol{  
  
  template <>
  void AutocorrelationStats<double>::append(const double& newValue, size_t iOfIteration, size_t iOfObservable, const double& newRefValue)
  {
    if (iOfIteration % decorrelationNbOfIterations_ == 0)
    {
      statisticsRefValues_.append(newRefValue, iOfObservable);
      indexRef_ = iOfIteration;
    }
    statisticsValues_.append(newValue, iOfObservable);
		//On intègre le profil d'autocorrélation à l'ordre 2
    statisticsMeanCorrelation_.append(((indexRef_ == iOfIteration)?.5:1) * statisticsRefValues_.lastValue(iOfObservable) * newValue, iOfObservable);
    statisticsCorrelation_.append(statisticsRefValues_.lastValue(iOfObservable) * newValue, ((iOfIteration - indexRef_)*nbOfAutocoPts_) / decorrelationNbOfIterations_, iOfObservable);    
  }
  
  template <>
  void AutocorrelationStats<dvec>::append(dvec const& newValue, size_t iOfIteration, size_t iOfObservable, dvec const& /*newRefValue*/)
  {
    if (iOfIteration % decorrelationNbOfIterations_ == 0)
    {
      statisticsRefValues_.append(newValue, iOfObservable);
      indexRef_ = iOfIteration;
    }
    statisticsValues_.append(newValue, iOfObservable);
		//On intègre le profil d'autocorrélation à l'ordre 2
    statisticsMeanCorrelation_.append(((indexRef_ == iOfIteration)?.5:1) * dot(statisticsRefValues_.lastValue(iOfObservable), newValue), iOfObservable);
    statisticsCorrelation_.append(dot(statisticsRefValues_.lastValue(iOfObservable), newValue), ((iOfIteration - indexRef_)*nbOfAutocoPts_) / decorrelationNbOfIterations_, iOfObservable);
  }
  
  template <>
  double AutocorrelationStats<double>::integratedAutocorrelationUnbiased(size_t iOfObservable) const
  {
    return max(0., integratedAutocorrelation(iOfObservable) - mean(iOfObservable) * statisticsRefValues_.mean(iOfObservable) * decorrelationTime_);
  }
  
  template <>
  double AutocorrelationStats<dvec>::integratedAutocorrelationUnbiased(size_t iOfObservable) const
  {
    return max(0., integratedAutocorrelation(iOfObservable) - dot(mean(iOfObservable), statisticsRefValues_.mean(iOfObservable)) * decorrelationTime_);
  }
  
}