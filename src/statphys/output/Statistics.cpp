#include "simol/statphys/output/Statistics.hpp"

using std::min;
using std::max;

namespace simol
{
  //##### CircularBuffer #####
  
  CircularBuffer::CircularBuffer(int size0):
    size_(size0),
    data_(DVec::Zero(size_)),
    index_(0),
    sum_(0),
    nbOfValues_(0)
    {}
    
    void CircularBuffer::append(double value)
    {
      if (nbOfValues_++ % (100 * size_) == 0)
      {
        //cout << "sum : " << sum() << " delta = " << sum() - data_.sum() << endl;
        data_(index_) = value;
        sum_ = data_.sum();
      }
      else
      {
        sum_ = sum_ + value - data_(index_);
        data_(index_) = value;
      }
      index_ = (index_+1)%size_;
      //if (value > 5)
      //  cout << data_.transpose() << endl;
    }
  
  //###### Statistics ######
  
  Statistics::Statistics(int nbRows, int nbCols):
    sumValues_(DMat::Zero(nbRows, nbCols)),
    lastValue_(DMat::Zero(nbRows, nbCols)),
    nbValues_(DMatInt::Zero(nbRows, nbCols))
    //iidVar_(nbRows, nbCols)
  {
    //sumValues_.fill(0);
    //lastValue_.fill(0);
    //nbValues_.fill(0);
    //iidVar_.fill(0);
  };
  
  int Statistics::nbOfRows() const
  {
    return sumValues_.rows();
  }
  
  int Statistics::nbOfCols() const
  {
    return sumValues_.cols();
  }

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
    sumVar_(DMat::Zero(nbRows, nbCols))
  {}
  
  void IIDStats::append(double value, int i, int j)
  {
    lastValue_(i, j) = value;
    nbValues_(i, j)++;
    double prevMean = mean();
    sumValues_(i, j) += value;
    sumVar_(i,j) = sumVar_(i,j) + (value - prevMean) * (value - mean());
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
    CorrelationStats(0,0,0,0)
    /*decorrelationNbOfSteps_(0),
    timeStep_(0),
    nbOfAutocoPts_(0),
    nbOfObservables_(0),
    statsValues_(0),
    statsRefValues_(0),
    statsIntegratedCorrelation_(0),
    statsCorrelation_(0, 0),
    indexRef_(0),
    bufferValues_(decorrelationNbOfSteps_)*/
  {
  }


  CorrelationStats::CorrelationStats(int decorrelationNbOfSteps0, double timeStep0, int nbOfAutocoPts0, int nbOfObservables):
    decorrelationNbOfSteps_(decorrelationNbOfSteps0),
    timeStep_(timeStep0),
    nbOfAutocoPts_(nbOfAutocoPts0),
    nbOfObservables_(nbOfObservables),
    statsValuesA_(nbOfObservables),
    statsValuesB_(nbOfObservables),
    statsIntegratedCorrelation_(nbOfObservables),
    statsCorrelation_(nbOfAutocoPts_, nbOfObservables),
    indexRef_(0),
    refValueB_(DVec::Zero(nbOfObservables)),
    currentCorrelationIntegral_(0),
    bufferValues_(decorrelationNbOfSteps_)
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
  
  void CorrelationStats::append(const double& newValueA, long int iOfStep, int iOfObservable, const double& newValueB)
  {
    //-- update the average --
    statsValuesA_.append(newValueA, iOfObservable);
    statsValuesB_.append(newValueB, iOfObservable);
    //-- update the autocorrelations --
    if (decorrelationNbOfSteps_ != 0)
    {
      if (iOfStep % decorrelationNbOfSteps_ == 0)
      {
        refValueB_(iOfObservable) = newValueB;
        indexRef_ = iOfStep;
        //cout << newValue << endl;
      }
      //-- update the correlation function --
      // decorrelationNbOfSteps_ / nbOfAutocoPts_ is the time intervale between two correlation estimations: can be bigger than the time step
      statsCorrelation_.append(refValueB_(iOfObservable) * newValueA, ((iOfStep - indexRef_)*nbOfAutocoPts_) / decorrelationNbOfSteps_, iOfObservable);
      
      //############################## OLD ESTIMATOR ##################
      /*//-- integration of the correlation function with a trapezoidal rule --
      // /!\ When the integral is small due the cancellation, this half factor can change the result up to 30% ! In this case dt is too large...
      //currentCorrelationIntegral_ += ((indexRef_ == iOfStep) ? .5 : 1) * statsRefValues_.lastValue(iOfObservable) * newValue;
      currentCorrelationIntegral_ += statsRefValues_.lastValue(iOfObservable) * newValue;
      if ((iOfStep+1) % decorrelationNbOfSteps_ == 0)    
      {
        statsIntegratedCorrelation_.append(timeStep() * currentCorrelationIntegral_, iOfObservable);
        currentCorrelationIntegral_ = 0;
      }*/
      
      //############################## NEW ESTIMATOR ##################
      bufferValues_.append(newValueA);
      if (bufferValues_.isReady())
        statsIntegratedCorrelation_.append(newValueB * 2 * bufferValues_.sumTrapeze() * timeStep());      // the factor 2 comes from the fact that we compute the integral on [-tdeco, tdeco]
        
      //cout << newValueB << " x " << bufferValues_.sumTrapeze() << endl;
    }
  }

  double CorrelationStats::correlationAtSpan(long int iOfSpan, int iOfObservable) const
  {
    return statsCorrelation_.mean(iOfSpan, iOfObservable);
  }
  
  double CorrelationStats::centeredCorrelationAtSpan(long int iOfSpan, int iOfObservable) const
  {
    return statsCorrelation_.mean(iOfSpan, iOfObservable) - meanA(iOfObservable) * statsValuesB_.mean(iOfObservable);
  }


  const double& CorrelationStats::lastValueA(int iOfObservable) const
  {
    return statsValuesA_.lastValue(iOfObservable);
  }


  double CorrelationStats::meanA(int iOfObservable) const
  {
    return statsValuesA_.mean(iOfObservable);
  }
  
  const long int& CorrelationStats::nbValues(int i, int j) const
  {
    return statsValuesA_.nbValues(i,j);
  }


  int CorrelationStats::statsIntegratedCorrelationNbValues(int iOfObservable) const
  {
    return statsIntegratedCorrelation_.nbValues(iOfObservable);
  }

  // This is the correlation integrated from -tdeco to tdeco
  double CorrelationStats::integratedCorrelation(int iOfObservable) const
  {
    //return statsIntegratedCorrelation_.mean(iOfObservable);
    return statsIntegratedCorrelation_.mean(iOfObservable);
  }

  // This is the correlation integrated from -tdeco to tdeco
  DVec CorrelationStats::integratedCorrelationVec() const
  {
    return statsIntegratedCorrelation_.meanMat().col(0);
  }

  // In the autocorrelation case this is the asymptotic variance
  double CorrelationStats::integratedCorrelationCentered(int iOfObservable) const
  {
    //cout << integratedCorrelation(iOfObservable) << " - " << 2 * meanA(iOfObservable) * statsValuesB_.mean(iOfObservable) * decorrelationTime() << endl;
    return integratedCorrelation(iOfObservable) - meanA(iOfObservable) * statsValuesB_.mean(iOfObservable) * (2*decorrelationTime() - timeStep());
  }
  
  // Returns the asymptotic variance of the estimator of the integrated correlation
  // See http://www.stat.unc.edu/faculty/cji/Sokal.pdf for a justification
  double CorrelationStats::asyvarIntegratedCorrelationCentered(int iOfObservable) const
  {
    return 4 * decorrelationTime() * pow(integratedCorrelationCentered(iOfObservable) , 2);
  }
  
  ///
  /// Returns an asymptotic variance: should be divided by the simulation time to obtain an error bar
  double CorrelationStats::varCorrelationAtSpan(int iOfSpan, int iOfObservable) const
  {
    if (iOfSpan < statsCorrelation_.nbOfRows())
      return statsCorrelation_.variance(iOfSpan, iOfObservable) * decorrelationTime();
    else
      return std::numeric_limits<double>::quiet_NaN();
      //return 0;
  }
  
  double CorrelationStats::stdDevCorrelationAtSpan(int iOfSpan, int iOfObservable) const
  {
    return sqrt(varCorrelationAtSpan(iOfSpan, iOfObservable));
  }
  
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
  }
  
  double AutocorrelationStats::centeredCorrelationAtSpan(long int iOfSpan, int iOfObservable) const
  {
    if (iOfSpan < statsCorrelation_.nbOfRows())
      return statsCorrelation_.mean(iOfSpan, iOfObservable) - pow(meanA(iOfObservable), 2);
    else
      return std::numeric_limits<double>::quiet_NaN();
      //return 0;
  }
  
  /*double AutocorrelationStats::integratedCorrelationCentered(int iOfObservable) const
  {
    cout << "AutocorrelationStats::integratedCorrelationCentered -> " << integratedCorrelation(iOfObservable) - 2 * pow(mean(iOfObservable),2) * decorrelationTime() << endl;
    return integratedCorrelation(iOfObservable) - 2 * pow(mean(iOfObservable),2) * decorrelationTime();
  }*/

  double AutocorrelationStats::asymptoticVariance(int iOfObservable) const
  {
    return integratedCorrelationCentered(iOfObservable);
  }

  /*double AutocorrelationStats::stdDev(int iOfObservable) const
  {
    return sqrt(variance(iOfObservable));
  }*/
  
  ///
  ///Returns the variance of the estimator of the variance
  double AutocorrelationStats::asyvarOfAsyvar(int iOfObservable) const
  {
    return asyvarIntegratedCorrelationCentered(iOfObservable);// * decorrelationTime();
  }
  
  /*///
  ///Returns the variance of the estimator of the variance
  double AutocorrelationStats::stdDevOfVar(int iOfObservable) const
  {
    return sqrt(varOfVar(iOfObservable));
  }*/
  
  /*///
  ///Returns the variance of the estimator of the variance
  double AutocorrelationStats::stdErrorOfVariance(int iOfObservable) const
  {
    return stdDevOfVariance(iOfObservable) / sqrt(statsCorrelation_.nbValues(iOfObservable));
  }*/
  
}
