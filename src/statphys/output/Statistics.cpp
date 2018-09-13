#include "simol/statphys/output/Statistics.hpp"

using std::min;
using std::max;

namespace simol
{
  /// Data structure with a circular index
  /// Keeps in memory the last "size_" given values and their sum  
  CircularBuffer::CircularBuffer(int size0):
    size_(size0),
    data_(DVec::Zero(size_)),
    index_(0),
    sum_(0),
    nbOfValues_(0)
    {}
    
    void CircularBuffer::append(double value)
    {
      // Every 100 new values, recompute the sum from scratch for numerical stability
      if (nbOfValues_++ % (100 * size_) == 0)
      {
        data_(index_) = value;
        sum_ = data_.sum();
      }
      else
      {
        sum_ = sum_ + value - data_(index_);
        data_(index_) = value;
      }
      index_ = (index_+1)%size_;
    }
  
  /// Statistics of nbRows X nbCols data samples, independently
  /// For each sample, computes the sum on the fly
  Statistics::Statistics(int nbRows, int nbCols):
    sumValues_(DMat::Zero(nbRows, nbCols)),
    lastValue_(DMat::Zero(nbRows, nbCols)),
    nbValues_(DMatInt::Zero(nbRows, nbCols))
  {
  };
  
  int Statistics::nbOfRows() const
  {
    return sumValues_.rows();
  }
  
  int Statistics::nbOfCols() const
  {
    return sumValues_.cols();
  }

  void Statistics::append(double value, int iOfRow, int iOfCol)
  {
    lastValue_(iOfRow, iOfCol) = value;
    sumValues_(iOfRow, iOfCol) += value;
    nbValues_(iOfRow, iOfCol)++;
  }


  double Statistics::mean(int iOfRow, int iOfCol) const
  {
    if (nbValues_(iOfRow, iOfCol) == 0)
      return 0;
    else
      return sumValues_(iOfRow, iOfCol) / nbValues(iOfRow, iOfCol);
  }

  /// Returns the vector of the means of the column i
  DVec Statistics::meanVec(int iOfCol) const
  {
    return sumValues_.col(iOfCol).array() / nbValues_.cast<double>().col(iOfCol).array();
  }

  /// Returns the matrix of the means
  DMat Statistics::meanMat() const
  {
    return sumValues_.array() / nbValues_.cast<double>().array();
  }


  const long int& Statistics::nbValues(int iOfRow, int iOfCol) const
  {
    return nbValues_(iOfRow, iOfCol);
  }

  const DVecInt Statistics::nbValuesVec(int iOfCol) const
  {
    return nbValues_.col(iOfCol);
  }

  const DMatInt& Statistics::nbValuesMat() const
  {
    return nbValues_;
  }


  const double& Statistics::lastValue(int iOfRow, int iOfCol) const
  {
    return lastValue_(iOfRow, iOfCol);
  }

  const DVec Statistics::lastValueVec(int iOfCol) const
  {
    return lastValue_.col(iOfCol);
  }

  const DMat& Statistics::lastValueMat() const
  {
    return lastValue_;
  }  
  
  
  /// Inherits from Statistics
  /// Computes the variance of each data sample on the fly, assuming they are IID  
  IIDStats::IIDStats(int nbRows, int nbCols):
    Statistics(nbRows, nbCols),
    sumVar_(DMat::Zero(nbRows, nbCols))
  {}
  
  void IIDStats::append(double value, int iOfRow, int iOfCol)
  {
    lastValue_(iOfRow, iOfCol) = value;
    nbValues_(iOfRow, iOfCol)++;
    double prevMean = mean();
    sumValues_(iOfRow, iOfCol) += value;
    sumVar_(iOfRow, iOfCol) = sumVar_(iOfRow, iOfCol) + (value - prevMean) * (value - mean());
  }
  
  double IIDStats::variance(int iOfRow, int iOfCol) const
  {
    if (nbValues_(iOfRow, iOfCol) == 1)
      return 0;
    else
      return sumVar_(iOfRow, iOfCol) / (nbValues_(iOfRow, iOfCol)-1.);
  }
  
  double IIDStats::stdDev(int iOfRow, int iOfCol) const
  {
    return sqrt(variance(iOfRow, iOfCol));
  }
  
  
  /// This class deals with two temporal series (A_n, B_n) of the form (phi(X_n), psi(X_n)) 
  /// Computes the correlation profile and its integrale
  CorrelationStats::CorrelationStats():
    CorrelationStats(0,0,0,0)
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
  }
  
  const int& CorrelationStats::decorrelationNbOfSteps() const
  {
    return decorrelationNbOfSteps_;
  }
  
  /// Expected upper bound on the time tau such that A_t and B_{t+tau} are independent
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
  
  const long int& CorrelationStats::nbValues(int iOfSpan, int iOfObservable) const
  {
    return statsValuesA_.nbValues(iOfSpan, iOfObservable);
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
  
  /// Specialization of CorrelationStats in the case when A = B
  //m Computes the autocorrelation profile and its integrale, which is the asymptotic variance
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

  double AutocorrelationStats::asymptoticVariance(int iOfObservable) const
  {
    return integratedCorrelationCentered(iOfObservable);
  }
  
  ///
  ///Returns the variance of the estimator of the variance
  double AutocorrelationStats::asyvarOfAsyvar(int iOfObservable) const
  {
    return asyvarIntegratedCorrelationCentered(iOfObservable);// * decorrelationTime();
  }
  
}
