#include "simol/statphys/output/Statistics.hpp"

using std::min;
using std::max;

namespace simol
{
  
  //###### Statistics ######
  
  Statistics::Statistics(int nbRows, int nbCols):
    sumValues_(nbRows, nbCols),
    lastValue_(nbRows, nbCols),
    nbValues_(nbRows, nbCols),
    iidVar_(nbRows, nbCols)
  {
    sumValues_.fill(0);
    lastValue_.fill(0);
    nbValues_.fill(0);
    iidVar_.fill(0);
  };

  void Statistics::append(double value, int i, int j)
  {
    lastValue_(i, j) = value;
    sumValues_(i, j) += value;
    nbValues_(i, j)++;
  }


  double Statistics::mean(int i, int j) const
  {
    if (nbValues_(i, j) == 0)
      return 0;
    else
      return sumValues_(i, j) / nbValues(i, j);
  }

  DVec Statistics::meanVec(int i) const
  {
    //return piecewiseDivision(sumValues_.col(i), nbValues_.col(i));
    return sumValues_.col(i).array() / nbValues_.cast<double>().col(i).array();
  }

  DMat Statistics::meanMat() const
  {
    //return piecewiseDivision(sumValues_, nbValues_);
    return sumValues_.array() / nbValues_.cast<double>().array();
  }


  const long int& Statistics::nbValues(int i, int j) const
  {
    return nbValues_(i, j);
  }

  const DVecInt Statistics::nbValuesVec(int i) const
  {
    return nbValues_.col(i);
  }

  const DMatInt& Statistics::nbValuesMat() const
  {
    return nbValues_;
  }


  const double& Statistics::lastValue(int i, int j) const
  {
    return lastValue_(i, j);
  }

  const DVec Statistics::lastValueVec(int i) const
  {
    return lastValue_.col(i);
  }

  const DMat& Statistics::lastValueMat() const
  {
    return lastValue_;
  }  
  
  
  //###### IIDStats ######
  
  IIDStats::IIDStats(int nbRows, int nbCols):
    Statistics(nbRows, nbCols),
    sumVar_(nbRows, nbCols)
  {
    sumVar_.fill(0);
  };
  
  void IIDStats::append(double value, int i, int j)
  {
    lastValue_(i, j) = value;
    nbValues_(i, j)++;
    double prevMean = mean();
    sumValues_(i, j) += value;
    sumVar_(i,j) = sumVar_(i,j) + (value - prevMean) * (value - mean());    
    
    /*lastValue_(i, j) = value;
    double tempTerm = iidVar_(i,j) + pow(mean(), 2);
    nbValues_(i, j)++;
    sumValues_(i, j) += value;
    iidVar_(i,j) = tempTerm - pow(mean(), 2) + (pow(value, 2) - tempTerm) / nbValues_(i, j);*/
  }
  
  double IIDStats::variance(int i, int j) const
  {
    if (nbValues_(i, j) == 1)
      return 0;
    else
    {
      //cout << sumVar_(i, j) << " / " << nbValues_(i, j)-1 << " = " << sumVar_(i, j) / (nbValues_(i, j)-1) << endl;
      return sumVar_(i, j) / (nbValues_(i, j)-1.);
    }
  }
  
  double IIDStats::stdDev(int i, int j) const
  {
    return sqrt(variance(i,j));
  }
  
  //###### CorrelationStats ######


  CorrelationStats::CorrelationStats():
    decorrelationNbOfSteps_(0),
    timeStep_(0),
    nbOfAutocoPts_(0),
    nbOfObservables_(0),
    //timeStep_(timeStep),
    statsValues_(0),
    statsRefValues_(0),
    statsIntegratedCorrelation_(0),
    //statsIntegratedCorrelationTemp_(),
    statsCorrelation_(0, 0),
    indexRef_(0)
  {
    //cout << "CorrelationStats()" << endl;
  }


  CorrelationStats::CorrelationStats(int decorrelationNbOfSteps0, double timeStep0, int nbOfAutocoPts0, int nbOfObservables):
    decorrelationNbOfSteps_(decorrelationNbOfSteps0),
    timeStep_(timeStep0),
    nbOfAutocoPts_(nbOfAutocoPts0),
    nbOfObservables_(nbOfObservables),
    //timeStep_(timeStep),
    statsValues_(nbOfObservables),
    statsRefValues_(nbOfObservables),
    statsIntegratedCorrelation_(nbOfObservables),
    //statsIntegratedCorrelationTemp_(),
    statsCorrelation_(nbOfAutocoPts_, nbOfObservables),
    indexRef_(0),
    currentCorrelationIntegral_(0)
  {
    //cout << "CorrelationStats(int decorrelation : nbObs = " << decorrelationTime << endl;
  }
  
  const int& CorrelationStats::decorrelationNbOfSteps() const
  {
    return decorrelationNbOfSteps_;
  }
  
  double CorrelationStats::decorrelationTime() const
  {
    return decorrelationNbOfSteps() * timeStep();
  }
  
  const double& CorrelationStats::timeStep() const
  {
    return timeStep_;
  }
  
  void CorrelationStats::append(const double& newValue, long int iOfStep, int iOfObservable, const double& newRefValue)
  {
    //-- update the average --
    statsValues_.append(newValue, iOfObservable);
    //-- update the autocorrelations --
    if (decorrelationNbOfSteps_ != 0)
    {
      if (iOfStep % decorrelationNbOfSteps_ == 0)
      {
        statsRefValues_.append(newRefValue, iOfObservable);
        indexRef_ = iOfStep;
        //cout << newValue << endl;
      }
      //-- update the correlation function --
      // decorrelationNbOfSteps_ / nbOfAutocoPts_ is the time intervale between two correlation estimations: can be bigger than the time step
      statsCorrelation_.append(statsRefValues_.lastValue(iOfObservable) * newValue, ((iOfStep - indexRef_)*nbOfAutocoPts_) / decorrelationNbOfSteps_, iOfObservable);
      //-- integration of the correlation function with a trapezoidal rule --
      // /!\ When the integral is small due the cancellation, this half factor can change the result up to 30% ! In this case dt is too big...
      //currentCorrelationIntegral_ += ((indexRef_ == iOfStep) ? .5 : 1) * statsRefValues_.lastValue(iOfObservable) * newValue;
      currentCorrelationIntegral_ += statsRefValues_.lastValue(iOfObservable) * newValue;
      //cout << currentCorrelationIntegral_ << endl;

      if ((iOfStep+1) % decorrelationNbOfSteps_ == 0)    
      {
        statsIntegratedCorrelation_.append(timeStep() * currentCorrelationIntegral_, iOfObservable);
        currentCorrelationIntegral_ = 0;
      }
    }
  }

  /*double CorrelationStats::operator()(long int iOfStep, int iOfObservable) const
  {
    return statsCorrelation_.mean(iOfStep, iOfObservable);
  }*/
  
  double CorrelationStats::correlationAtSpan(long int iOfSpan, int iOfObservable) const
  {
    return statsCorrelation_.mean(iOfSpan, iOfObservable);
  }
  
  double CorrelationStats::unbiasedCorrelationAtSpan(long int iOfSpan, int iOfObservable) const
  {
    return statsCorrelation_.mean(iOfSpan, iOfObservable) - mean(iOfObservable) * statsRefValues_.mean(iOfObservable);
  }


  const double& CorrelationStats::lastValue(int iOfObservable) const
  {
    return statsValues_.lastValue(iOfObservable);
  }


  double CorrelationStats::mean(int iOfObservable) const
  {
    return statsValues_.mean(iOfObservable);
  }
  
  const long int& CorrelationStats::nbValues(int i, int j) const
  {
    return statsValues_.nbValues(i,j);
  }


  int CorrelationStats::statsIntegratedCorrelationNbValues(int iOfObservable) const
  {
    return statsIntegratedCorrelation_.nbValues(iOfObservable);
  }


  double CorrelationStats::integratedCorrelation(int iOfObservable) const
  {
    return statsIntegratedCorrelation_.mean(iOfObservable);
  }

  DVec CorrelationStats::integratedCorrelationVec() const
  {
    return statsIntegratedCorrelation_.meanMat().col(0);
  }

  double CorrelationStats::integratedCorrelationUnbiased(int iOfObservable) const
  {
    return integratedCorrelation(iOfObservable) - mean(iOfObservable) * statsRefValues_.mean(iOfObservable) * decorrelationTime();
  }
  
  double CorrelationStats::varIntegratedCorrelation(int iOfObservable) const
  {
    return statsIntegratedCorrelation_.variance(iOfObservable);
  }
  
  ///
  /// Returns an asymptotic variance: should be divided by the simulation time to obtain an error bar
  double CorrelationStats::varCorrelationAtSpan(int iOfSpan, int iOfObservable) const
  {
    return statsCorrelation_.variance(iOfSpan, iOfObservable) * decorrelationTime();
  }
  
  double CorrelationStats::stdDevCorrelationAtSpan(int iOfSpan, int iOfObservable) const
  {
    return sqrt(varCorrelationAtSpan(iOfSpan, iOfObservable));
  }
  
  /*double CorrelationStats::stdErrorCorrelationAtSpan(int iOfSpan, int iOfObservable) const
  {
    return stdDevCorrelationAtSpan(iOfSpan, iOfObservable) / sqrt(statsCorrelation_.nbValues(iOfObservable));
  }*/
  
  //###### AutocorrelationStats ######
  
  AutocorrelationStats::AutocorrelationStats():
    CorrelationStats()
  {}


  AutocorrelationStats::AutocorrelationStats(int decorrelationNbOfSteps0, double timeStep0, int nbOfAutocoPts0, int nbOfObservables):
    CorrelationStats(decorrelationNbOfSteps0, timeStep0, nbOfAutocoPts0, nbOfObservables)
  {}
  
  void AutocorrelationStats::append(double const& newValue, long int iOfStep, int iOfObservable)
  {
    CorrelationStats::append(newValue, iOfStep, iOfObservable, newValue);
      if ((iOfStep+1) % decorrelationNbOfSteps_ == 0)    
      {        
        /*cout << "c0 = " << correlationAtSpan(0)<< endl;
        double test = 0;
        for (int i=0; i < nbOfAutocoPts_; i++)
        {
          test += (i==0?.5:1)*unbiasedCorrelationAtSpan(i) * timeStep() * decorrelationNbOfSteps() / nbOfAutocoPts_;
          //cout << "+= " << 2 * unbiasedCorrelationAtSpan(i) * timeStep() * decorrelationNbOfSteps() / nbOfAutocoPts_ << endl;
        }
        cout << "---->" << integratedCorrelationUnbiased() << " <-> " << test - .5  * unbiasedCorrelationAtSpan(0)* timeStep() * decorrelationNbOfSteps() / nbOfAutocoPts_ << " = " << test << " - " << .5 * unbiasedCorrelationAtSpan(0)* timeStep() * decorrelationNbOfSteps() / nbOfAutocoPts_ <<  endl;
        //cout << "---->" << variance() << " <-> " << test - .5 * 2 * unbiasedCorrelationAtSpan(0)* timeStep() * decorrelationNbOfSteps() / nbOfAutocoPts_ << " = " << test << " - " << .5 * 2 * unbiasedCorrelationAtSpan(0, 0)* timeStep() * decorrelationNbOfSteps() / nbOfAutocoPts_ <<  endl;
      */
      }
  }
  
  double AutocorrelationStats::unbiasedCorrelationAtSpan(long int iOfSpan, int iOfObservable) const
  {
    return statsCorrelation_.mean(iOfSpan, iOfObservable) - pow(mean(iOfObservable), 2);
  }
  
  double AutocorrelationStats::integratedCorrelationUnbiased(int iOfObservable) const
  {
    return integratedCorrelation(iOfObservable) - pow(mean(iOfObservable),2) * decorrelationTime();
  }

  double AutocorrelationStats::variance(int iOfObservable) const
  {
    //return 2 * max(0., integratedCorrelationUnbiased(iOfObservable));
    return 2 * integratedCorrelationUnbiased(iOfObservable);
  }

  double AutocorrelationStats::stdDev(int iOfObservable) const
  {
    return sqrt(variance(iOfObservable));
  }
  
  ///
  ///Returns the variance of the estimator of the variance
  double AutocorrelationStats::varOfVar(int iOfObservable) const
  {
    return 4 * varIntegratedCorrelation(iOfObservable) * decorrelationTime();
  }
  
  ///
  ///Returns the variance of the estimator of the variance
  double AutocorrelationStats::stdDevOfVar(int iOfObservable) const
  {
    return sqrt(varOfVar(iOfObservable));
  }
  
  /*///
  ///Returns the variance of the estimator of the variance
  double AutocorrelationStats::stdErrorOfVariance(int iOfObservable) const
  {
    return stdDevOfVariance(iOfObservable) / sqrt(statsCorrelation_.nbValues(iOfObservable));
  }*/
  
}
