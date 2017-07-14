#include "simol/statphys/controlVariate/ControlVariate.hpp"

using std::cout;
using std::endl;
using std::sin;
using std::cos;
using std::exp;
using std::ofstream;

namespace simol
{

  Observable* createControlVariate(Input const& input, int idObs, shared_ptr<CVBasis> cvBasis0)
  {
    if (input.controlVariateName() == "None")
      return new Observable(input, idObs, input.decorrelationNbOfSteps(), input.nbOfAutocoPts());
    if (input.controlVariateName() == "Sinus")
      return new SinusControlVariate(input, idObs, cvBasis0);
    else if (input.controlVariateName() == "Basis")
      return new BasisControlVariate(input, idObs, cvBasis0);
    else
      std::cout << input.controlVariateName() << " is not a valid control variate !" << std::endl;
    return 0; 
    
      /*ExpFourierHermiteControlVariate(input, idObs, cvBasis0);
    else if(input.controlVariateName() == "HermiteHermite")
      return new HermiteHermiteControlVariate(input, idObs, cvBasis0);
    else if(input.controlVariateName() == "ExpHermiteHermite")
      return new ExpHermiteHermiteControlVariate(input, idObs, cvBasis0);
    else if(input.controlVariateName() == "QuadraticHermite")
      return new QuadraticHermiteControlVariate(input, idObs, cvBasis0);*/

  }


  ControlVariate::ControlVariate(Input const& input, int idObs, shared_ptr<CVBasis> cvBasis0, int nbOfFunctions):
    Observable(input, idObs, input.decorrelationNbOfSteps(), input.nbOfAutocoPts()),
    dimension_(input.dimension()),
    nbOfFunctions_(nbOfFunctions),
    nbOfFunctionPairs_(pow(nbOfFunctions_, 2)),
    doEstimateCvCoeffs_(false),
    CVDecorrelationNbOfSteps_(input.shortDecorrelationNbOfSteps()),
    CVNbOfAutocoPts_(input.nbOfShortAutocoPts()),
    autocoStatsBetter_(CVDecorrelationNbOfSteps(), timeStep(), CVNbOfAutocoPts()),
    statsGeneratorOnBasis_(nbOfFunctions_),
    statsB1_(nbOfFunctions_),
    statsB2_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfFunctions_),
    statsD_(nbOfFunctions_, nbOfFunctions_),
    lastA_(nbOfFunctions_),
    cvBasis_(cvBasis0)
  {    
    cout << "decorrelationNbOfSteps : " << input.decorrelationNbOfSteps() << endl;
    cout << "shortDecorrelationNbOfSteps : " << input.shortDecorrelationNbOfSteps() << endl;
    /*if (!cvBasis_) throw runtime_error("ControlVariate build without cvBasis !");
    if (cvBasis().cvCoeffs_)
    {
      cout << "CVcoeffs from Galerkin !" << endl;    
    }
    else if (input.controlVariateCoeffsPath() != "None")
    {
      cout << "CVcoeffs from file !" << endl;
      //throw runtime_error("Reading CV coefficients from a file is not fixed yet !");
      std::string coeffsPath = input.outputFolderName() + input.controlVariateCoeffsPath();
      //ifstream file(coeffsPath);
      vector<int> dimensions;
      cvBasis().cvCoeffs_ = make_shared<DVec>(scanTensor(coeffsPath, dimensions));
      //coeffsVec_ = SMat(coeffsPath, getNbOfLines(file));
    }
    else
    {
      cout << "CVcoeffs estimated on the fly !" << endl;
      doEstimateCvCoeffs_ = true;
    }*/
    if (!cvBasis().cvCoeffs_)
      doEstimateCvCoeffs_ = true;
  }

  int ControlVariate::nbOfFunctions() const
  {
    return nbOfFunctions_;
  }

  int ControlVariate::nbOfFunctionPairs() const
  {
    return nbOfFunctionPairs_;
  }

  /*bool ControlVariate::doOutput(long int iOfStep) const
  {
    return (printPeriodNbOfSteps_ > 0 && iOfStep % printPeriodNbOfSteps_ == 0);
  }*/

  bool ControlVariate::isNone() const
  {
    return false;
  }
  
  bool ControlVariate::doEstimateCvCoeffs() const
  {return doEstimateCvCoeffs_;}


  int ControlVariate::nbOfFourier() const
  {return 0;}

  int ControlVariate::nbOfHermite() const
  {return 0;}
  
  int ControlVariate::CVDecorrelationNbOfSteps() const
  {return CVDecorrelationNbOfSteps_;}
  
  int ControlVariate::CVNbOfAutocoPts() const
  {return CVNbOfAutocoPts_;}
  
  
  
  

  double ControlVariate::lastValueBetter() const
  {
    return autocoStatsBetter_.lastValue();
  }

  double ControlVariate::meanBetter() const
  {
    return autocoStatsBetter_.mean();
  }
  
  double ControlVariate::varBetter() const
  {
    return autocoStatsBetter_.variance();
  }

  double ControlVariate::stdDevBetter() const
  {
    return autocoStatsBetter_.stdDev();
  }
  
  double ControlVariate::varOfVarBetter() const
  {
    return autocoStatsBetter_.varOfVar();
  }

  DVec ControlVariate::correlationB2() const
  {
    return statsB2_.integratedCorrelationVec();
  }

  double ControlVariate::correlationB2(int iOfFunction) const
  {
    return statsB2_.integratedCorrelation(iOfFunction);
  }
  
  double ControlVariate::correlationBetterAtSpan(int iOfSpan) const
  {
    return autocoStatsBetter_.correlationAtSpan(iOfSpan);
  }
  
  double ControlVariate::unbiasedCorrelationBetterAtSpan(int iOfSpan) const
  {
    //return autocoStats_(iOfSpan) - pow(meanObservable(), 2);
    return autocoStatsBetter_.unbiasedCorrelationAtSpan(iOfSpan);
  }
  
  double ControlVariate::varCorrelationBetterAtSpan(int iOfSpan) const
  {
    return autocoStatsBetter_.varCorrelationAtSpan(iOfSpan);
  }

  double ControlVariate::stdDevCorrelationBetterAtSpan(int iOfSpan) const
  {
    return autocoStatsBetter_.stdDevCorrelationAtSpan(iOfSpan);
  }
  
  
  
  
  CVBasis const& ControlVariate::cvBasis() const
  {return *cvBasis_;}
  
  CVBasis& ControlVariate::cvBasis()
  {return *cvBasis_;}
  
  const double& ControlVariate::basisValue(int iOfFunction) const
  {
    return cvBasis().basisValues_(iOfFunction);
  }
  
  const DVec& ControlVariate::basisValues() const
  {
    return cvBasis().basisValues_;
  }
  
  const double& ControlVariate::generatorOnBasisValue(int iOfFunction) const
  {
    return cvBasis().generatorOnBasisValues_(iOfFunction);
  }
  
  const DVec& ControlVariate::generatorOnBasisValues() const
  {
    return cvBasis().generatorOnBasisValues_;
  }
  
  

  void ControlVariate::appendToObservable(double observable, long int iOfStep)
  {
    autocoStats_.append(observable, iOfStep);
  }

  void ControlVariate::appendToB1(double observable)
  {
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
      statsB1_.append((observable - autocoStats_.mean()) * basisValue(iOfFunction), iOfFunction);
  }

  void ControlVariate::appendToB2(double observable, long int iOfStep)
  {
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
      statsB2_.append((observable - autocoStats_.mean()), iOfStep, iOfFunction, generatorOnBasisValue(iOfFunction));
  }

  void ControlVariate::appendToD()
  {
    //Eigen::Matrix<CorrelationStats<double>, Eigen::Dynamic, Eigen::Dynamic> A(nbOfFunctions_, nbOfFunctions_, CorrelationStats<double>(decorrelationNbOfSteps(), decorrelationTime()));

    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
      for (int iOfFunction2 = 0; iOfFunction2 <= iOfFunction; iOfFunction2++)
      {
        double valueSym = (- basisValue(iOfFunction) * generatorOnBasisValue(iOfFunction2)
                           - basisValue(iOfFunction2) * generatorOnBasisValue(iOfFunction)) / 2.;
        statsD_.append(valueSym, iOfFunction, iOfFunction2);
        if (iOfFunction != iOfFunction2)
          statsD_.append(valueSym, iOfFunction2, iOfFunction);
      }
  }



  void ControlVariate::appendToBetter(double observable, long int iOfStep)
  {    
    if (doEstimateCvCoeffs())
    {
      throw std::runtime_error("CV with coefficients estimated on the fly is not fixed yet !");
      //double betterObservableTerm = dot((statsB_.meanMat() - statsB2_.integratedCorrelationMat()), statsD_.meanMat().llt().solve(generatorOnBasisFunction));
      /*DMat Dinv = statsD_.meanMat().llt().solve(generatorOnBasisFunction);
      DVec B = statsB_.meanMat() - statsB2_.integratedCorrelationMat();
      std::cout << B.transpose() * Dinv << endl;*/
      //cout << (statsB_.meanMat() - statsB2_.integratedCorrelationMat()).transpose().size() << "  " << statsD_.meanMat().inverse().size() << endl;
      //lastA_ = (statsB_.meanMat() - statsB2_.integratedCorrelationMat()).transpose() * statsD_.meanMat().llt().solve(generatorOnBasisFunction);
      //cout << lastA_.size() << "   " << generatorOnBasisFunction.size() << endl;
    }
    else
    {
      if (!cvBasis().cvCoeffs_) throw std::runtime_error("cvCoeffs is supposed to be initialized here !");
      //cout << "---------> " << observable << " + " << dot(*cvBasis().cvCoeffs_, generatorOnBasisValues()) << endl;
      //cout << "---------> " << *cvBasis().cvCoeffs_ << " X " << generatorOnBasisValues() << endl;
      //cout << "ControlVariate::appendToBetter : " << observable + dot(*cvBasis().cvCoeffs_, generatorOnBasisValues()) << endl;
      autocoStatsBetter_.append(observable + dot(*cvBasis().cvCoeffs_, generatorOnBasisValues()), iOfStep);
    }

    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
    {
      //historyGeneratorOnBasis_(iOfStep, iOfFunction) = generatorOnBasisFunction(iOfFunction);
      statsGeneratorOnBasis_.append(generatorOnBasisValue(iOfFunction), iOfFunction);
    }
  }



  DVec ControlVariate::lastB1() const
  {
    return statsB1_.lastValueVec();
  }


  DVec ControlVariate::meanB1() const
  {
    return statsB1_.meanVec();
  }

  DVec ControlVariate::meanB() const
  {
    return statsB1_.meanVec() - statsB2_.integratedCorrelationVec();
  }

  DMat ControlVariate::lastD() const
  {
    return statsD_.lastValueMat();
  }

  DMat ControlVariate::meanD() const
  {
    return statsD_.meanMat();
  }

  DVec ControlVariate::lastGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.lastValueVec();
  }

  DVec ControlVariate::meanGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.meanVec();
  }

  double ControlVariate::autocorrelationB2(int iOfSpan, int iOfFunction) const
  {
    return statsB2_.correlationAtSpan(iOfSpan, iOfFunction);
  }

  DVec ControlVariate::lastA() const
  {
    return lastA_;

  }
  double ControlVariate::lastA(int iOfFunction) const
  {
    return lastA_(iOfFunction);
  }

  ///
  ///Feed a new value to the Observable
  /// /!\ an update*** method must be called beforehand
  void ControlVariate::append(double observable, long int iOfStep)
  {
    if (doEstimateCvCoeffs())
    {
      appendToB1(observable);
      appendToB2(observable, iOfStep);
      appendToD();
    }
    appendToObservable(observable, iOfStep);
    appendToBetter(observable, iOfStep);
    //cout << "end ControlVariate::update" << endl;
  }
  


  void ControlVariate::display(long int iOfStep)
  {
    outFlux() << iOfStep * timeStep();
    if (doEstimateCvCoeffs())
    {
      for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
        outFlux() << " " << meanB()(iOfFunction) * lastA(iOfFunction);
      for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
        outFlux() << " " << lastA(iOfFunction);
      /*for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
        outFlux() << " " << lastB()(iOfFunction);*/
      for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
        outFlux() << " " << meanB()(iOfFunction);
      /*for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
        for (int iOfFunction2 = 0; iOfFunction2 < nbOfFunctions(); iOfFunction2++)
      outFlux() << " " << lastD()(iOfFunction, iOfFunction2);*/
      for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
        for (int iOfFunction2 = 0; iOfFunction2 < nbOfFunctions(); iOfFunction2++)
          outFlux() << " " << meanD()(iOfFunction, iOfFunction2);  //8-11
    }
    outFlux() << " " << lastValue()
        << " " << mean()   //13
        << " " << variance()
        << " " << varOfVar()
        << " " << lastValueBetter()
        << " " << meanBetter()
        << " " << varBetter()
        << " " << varOfVarBetter()
        << " " << dot(*cvBasis().cvCoeffs_, basisValues());
    if (doEstimateCvCoeffs())
    {
      for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
        outFlux() << " " << lastGeneratorOnBasis()(iOfFunction);
      for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
        outFlux() << " " << meanGeneratorOnBasis()(iOfFunction);

      for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
        outFlux() << " " << meanB1()(iOfFunction);

      for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
        outFlux() << " " << -correlationB2()(iOfFunction);  //22-23
    }
    outFlux() << endl;
  }
  
  void ControlVariate::displayFinalValues(ofstream& out)
  {
    out << mean()
        << " " << variance()
        << " " << varOfVar() 
        << " " << meanBetter()
        << " " << varBetter()
        << " " << varOfVarBetter()
        << endl;
  }
  
  void ControlVariate::displayCorrelations(long int iOfStep)
  {
    //cout << "Velocity : The correlation in 0 is " << floor((2 * obsVelocity().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsVelocity().variance() * 10000)/100 << "% of the variance" << endl;
    //cout << "SumFlow : The correlation in 0 is " << floor((2 * obsSumFlow().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsSumFlow().variance() * 10000)/100 << "% of the variance" << endl;
    //cout << "ModiFlow : The correlation in 0 is " << floor((2 * obsModiFlow().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsModiFlow().variance() * 10000)/100 << "% of the variance" << endl;
       
    //cout << "Velocity : The correlation in 0 is " << floor((2 * obsVelocity().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsVelocity().variance() * 10000)/100 << "% of the variance" << endl;
    outFluxCorrelation() << "# The correlation in iOfSpan = 0 is " << floor((2 * unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / variance() * 10000)/100 << "% of the variance" << endl;
    outFluxCorrelation() << "# The correlation in iOfSpan = 0 is " << floor((2 * unbiasedCorrelationBetterAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / varBetter() * 10000)/100 << "% of the better variance" << endl;

    for (int iOfSpan = 0; iOfSpan < nbOfAutocoPts(); iOfSpan++)
    {
      outFluxCorrelation() << iOfSpan * autocoPtsPeriod()
                       << " " << unbiasedCorrelationAtSpan(iOfSpan)
                       << " " << sqrt(varCorrelationAtSpan(iOfSpan) / (iOfStep * timeStep()))
                       << " " << unbiasedCorrelationBetterAtSpan(iOfSpan)
                       << " " << sqrt(varCorrelationBetterAtSpan(iOfSpan) / (iOfStep * timeStep())) << endl;
    }
    
    /*for (int iOfSpan = nbOfAutocoPts(); iOfSpan < nbOfAutocoPts(); iOfSpan++)
    {
      outFluxCorrelation() << iOfSpan * autocoPtsPeriod()
                       << " " << unbiasedCorrelationAtSpan(iOfSpan)
                       << " " << sqrt(varCorrelationAtSpan(iOfSpan) / (iOfStep * timeStep()))
                       << " " << 0
                       << " " << 0 << endl;
    }*/
  }



 

  SinusControlVariate::SinusControlVariate(Input const& input, int idObs, shared_ptr<CVBasis> cvBasis0):
    ControlVariate(input, idObs, cvBasis0, 1)
  {}

  double SinusControlVariate::value(System const& syst) const
  {
    double q = syst(0).position(0);
    return sin(2 * M_PI * q);
  }

  /*DVec SinusControlVariate::gradientQ(System const& syst, int iOfParticle, int iOfFunction) const
  {
    double q = syst(0).position(0);
    DVec grad(dimension_);
    grad(0) = 2 * M_PI * cos(2 * M_PI * q);
    return grad;
  }

  double SinusControlVariate::laplacianQ(System const& syst, int iOfParticle, int iOfFunction) const
  {
    double q = syst(0).position(0);
    return - pow (2 * M_PI, 2) * sin(2 * M_PI * q);
  }


  double SinusControlVariate::laplacianP(System const&, int iOfParticle, int iOfFunction) const
  {
    return 0;
  }

  DVec SinusControlVariate::gradientP(System const&, int iOfParticle, int iOfFunction) const
  {
    return DVec(dimension_);
  }*/


  //#####BasisControlVariate#####

  BasisControlVariate::BasisControlVariate(Input const& input, int idObs, shared_ptr<CVBasis> cvBasis0):
    ControlVariate(input, idObs, cvBasis0, 1)
  {

    //cout << "coeffsVec_ : " << endl;
    //cout << coeffsVec_ << endl;
  }
  
  DVec const& BasisControlVariate::cvCoeffs() const
  {
    if (!(cvBasis().cvCoeffs_)) throw runtime_error("cvCoeffs not in memory !");
    return (*cvBasis().cvCoeffs_);
  }
  
  DVec& BasisControlVariate::cvCoeffs()
  {
    if (!(cvBasis().cvCoeffs_)) throw runtime_error("cvCoeffs not in memory !");
    return (*cvBasis().cvCoeffs_);
  }
  
  double const& BasisControlVariate::cvCoeffs(int i) const
  {
    if (!(cvBasis().cvCoeffs_)) throw runtime_error("cvCoeffs not in memory !");
    return (*cvBasis().cvCoeffs_)(i);
  }
  
  double& BasisControlVariate::cvCoeffs(int i)
  {
    if (!(cvBasis().cvCoeffs_)) throw runtime_error("cvCoeffs not in memory !");
    return (*cvBasis().cvCoeffs_)(i);
  }

  double BasisControlVariate::value(System const&) const
  {
    //cout << "BasisControlVariate::value : " << dot(cvCoeffs(), basisValues()) << endl;
    return dot(cvCoeffs(), basisValues());
  }
  
   void BasisControlVariate::display(long int iOfStep)
  {
    outFlux() << iOfStep * timeStep()
        << " " << lastValue()
        << " " << mean()
        << " " << variance()
        << " " << varOfVar()
        << " " << lastValueBetter()
        << " " << meanBetter()
        << " " << varBetter()
        << " " << varOfVarBetter()
        << " " << dot(*cvBasis().cvCoeffs_, basisValues())
        << endl;
  }

  /*DMat BasisControlVariate::gradientQ(System const& syst, int iOfParticle, int iOfFunction) const
  {
    //cout << "BasisControlVariate::gradientQ" << endl;
    DMat result = DMat::Zero(syst.dimension(), syst.nbOfParticles());
    
    for (int iOfFunction=0; iOfFunction < (int)cvCoeffs().size(); iOfFunction++)
      result += cvCoeffs(iOfFunction) * cvBasis().basis_->gradientQ(syst, iOfFunction);
      
    return result;
  }

  double BasisControlVariate::laplacianQ(System const& syst, int iOfFunction) const
  {
    double result = 0;

    for (int iOfFunction=0; iOfFunction < (int)cvCoeffs().size(); iOfFunction++)
      result += cvCoeffs(iOfFunction) * cvBasis().laplacianQ(syst, iOfFunction);

    return result;
  }

  DVec BasisControlVariate::gradientP(System const& syst, int iOfFunction) const
  {
    DVec result(1, 0);
    
    for (int i=0; i < (int)cvCoeffs().size(); i++)
      result += cvCoeffs(i) * cvBasis().basis_->gradientP(syst, iOfParticle, i);
    //cout << endl;
    //cout << "gradientP = " << result << endl;
    return result;
  }

  double BasisControlVariate::laplacianP(System const& syst, int iOfParticle, int iOfFunction) const
  {
    double result = 0;

    for (int i=0; i < (int)cvCoeffs().size(); i++)
      result += cvCoeffs(i) * cvBasis().basis_->laplacianP(syst, iOfParticle, i);
    //cout << "laplacianP = " << result << endl;
    return result;
  }*/



  ExpFourierHermiteControlVariate::ExpFourierHermiteControlVariate(Input const& input, int idObs, shared_ptr<CVBasis> cvBasis0):
    BasisControlVariate(input, idObs, cvBasis0)
  {
    //basis_ = new ExpFourierHermiteBasis(input, potential);
    nbQ_ = 100;  // rafinement du maillage de l'output map
    nbP_ = 100;  // rafinement du maillage de l'output map
    deltaQ_ = 2 * M_PI / nbQ_;
    pMax_ = 4;
    deltaP_ = 2 * pMax_ / nbP_;
  }

  int ExpFourierHermiteControlVariate::nbOfFourier() const
  {return cvBasis().tensorBasis_->nbOfElts(0);}

  int ExpFourierHermiteControlVariate::nbOfHermite() const
  {return cvBasis().tensorBasis_->nbOfElts(1);}

  /*void ExpFourierHermiteControlVariate::displayMap(ofstream& out) const
  {
    cout << "displayMap" << endl;
    vector<Particle*> conf(1, new Particle(0, 0, 0));
    //conf[0] = Particle(0,0,0);
    out << nbP_ + 1 << " ";
    for (double p = -pMax_; p < pMax_; p += deltaP_)
      out << p << " ";
    out << endl;
    for (double q = -M_PI; q < M_PI; q += deltaQ_)
    {
      out << q << " ";
      for (double p = -pMax_; p < pMax_; p += deltaP_)
      {
        conf[0]->position() = DVec::Constant(1, q);
        conf[0]->momentum() = DVec::Constant(1, p);
        out << basisFunction(conf, 0) << " ";
      }
      out << endl;
    }
    cout << "end displayMap" << endl;
  }*/
  
  /*HermiteHermiteControlVariate::HermiteHermiteControlVariate(Input const& input, int idObs, shared_ptr<CVBasis> cvBasis0):
    BasisControlVariate(input, idObs, cvBasis0)
  {}
  
  ExpHermiteHermiteControlVariate::ExpHermiteHermiteControlVariate(Input const& input, int idObs, shared_ptr<CVBasis> cvBasis0):
    BasisControlVariate(input, idObs, cvBasis0)
  {}
  
  QuadraticHermiteControlVariate::QuadraticHermiteControlVariate(Input const& input, int idObs, shared_ptr<CVBasis> cvBasis0):
    BasisControlVariate(input, idObs, cvBasis0)
  {}*/

}










