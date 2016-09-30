#include "simol/statphys/controlVariate/ControlVariate.hpp"

using std::cout;
using std::endl;
using std::sin;
using std::cos;
using std::exp;
using std::ofstream;

namespace simol
{

  Observable* createControlVariate(Input const& input, const string& outPath, CVBasis& cvBasis0, Galerkin* galerkin)
  {
    if (input.controlVariateName() == "None")
      return new Observable(input, outPath);
    if (input.controlVariateName() == "Sinus")
      return new SinusControlVariate(input, outPath, cvBasis0);
    else if (input.controlVariateName() == "ExpFourierHermite")
      return new ExpFourierHermiteControlVariate(input, outPath, cvBasis0, galerkin);
    else
      std::cout << input.controlVariateName() << " is not a valid control variate !" << std::endl;
    return 0;
  }


  ControlVariate::ControlVariate(Input const& input, const string& outPath, CVBasis& cvBasis0, int nbOfFunctions):
    Observable(input, outPath),
    dimension_(input.dimension()),
    nbOfFunctions_(nbOfFunctions),
    nbOfFunctionPairs_(pow(nbOfFunctions_, 2)),
    autocoStatsBetter_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts()),
    statsGeneratorOnBasis_(nbOfFunctions_),
    statsB1_(nbOfFunctions_),
    statsB2_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfFunctions_),
    statsD_(nbOfFunctions_, nbOfFunctions_),
    lastA_(nbOfFunctions_),
    cvBasis_(&cvBasis0)
    /*,
    basisValue(nbOfFunctions_),
    generatorOnBasisValue(nbOfFunctions_)*/
  {}

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

  /*double ControlVariate::potential(Vector<double> const& position) const
  {
    return (*potential_)(position);
  }

  Vector<double> ControlVariate::potentialDerivative(Vector<double> const& position) const
  {
    return potential_->gradient(position);
  }

  double ControlVariate::potentialLaplacian(Vector<double> const& position) const
  {
    return potential_->laplacian(position);
  }*/

  int ControlVariate::nbOfFourier() const
  {return 0;}

  int ControlVariate::nbOfHermite() const
  {return 0;}

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

  Vector<double> ControlVariate::correlationB2() const
  {
    return statsB2_.integratedCorrelationVec();
  }

  double ControlVariate::correlationB2(int iOfFunction) const
  {
    return statsB2_.integratedCorrelation(iOfFunction);
  }
  
  const double& ControlVariate::basisValue(int iOfFunction) const
  {
    return cvBasis_->basisValues_(iOfFunction);
  }
  
  const DVec& ControlVariate::basisValues() const
  {
    return cvBasis_->basisValues_;
  }
  
  const double& ControlVariate::generatorOnBasisValue(int iOfFunction) const
  {
    return cvBasis_->generatorOnBasisValues_(iOfFunction);
  }
  
  const DVec& ControlVariate::generatorOnBasisValues() const
  {
    return cvBasis_->generatorOnBasisValues_;
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
    //double betterObservableTerm = dot((statsB_.meanMat() - statsB2_.integratedCorrelationMat()), statsD_.meanMat().llt().solve(generatorOnBasisFunction));
    /*DenseMatrix<double> Dinv = statsD_.meanMat().llt().solve(generatorOnBasisFunction);
    Vector<double> B = statsB_.meanMat() - statsB2_.integratedCorrelationMat();
    std::cout << B.transpose() * Dinv << endl;*/
    //cout << (statsB_.meanMat() - statsB2_.integratedCorrelationMat()).transpose().size() << "  " << statsD_.meanMat().inverse().size() << endl;
    //lastA_ = (statsB_.meanMat() - statsB2_.integratedCorrelationMat()).transpose() * statsD_.meanMat().llt().solve(generatorOnBasisFunction);
    //cout << lastA_.size() << "   " << generatorOnBasisFunction.size() << endl;

    if (statsD_.meanMat().determinant() != 0)
    {
      //lastA_ = - .5 * statsD_.meanMat().llt().solve(meanB());
      //lastA_ = Vector<double>(1,1);
      lastA_.fill(1);
      autocoStatsBetter_.append(observable - dot(lastA_, generatorOnBasisValues()), iOfStep);
    }
    else
    {
      lastA_.fill(0);
      autocoStatsBetter_.append(observable, iOfStep);
    }

    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
    {
      //historyGeneratorOnBasis_(iOfStep, iOfFunction) = generatorOnBasisFunction(iOfFunction);
      statsGeneratorOnBasis_.append(generatorOnBasisValue(iOfFunction), iOfFunction);
    }
  }



  Vector<double> ControlVariate::lastB1() const
  {
    /*Vector<double> result(nbOfFunctions_);
    for (int iOfFunction =0; iOfFunction < nbOfFunctions_; iOfFunction++)
      result(iOfFunction) = statsB_.lastValue(iOfFunction);
    return result;*/
    return statsB1_.lastValueVec();
  }


  Vector<double> ControlVariate::meanB1() const
  {
    return statsB1_.meanVec();
  }

  Vector<double> ControlVariate::meanB() const
  {
    return statsB1_.meanVec() - statsB2_.integratedCorrelationVec();
  }

  DenseMatrix<double> ControlVariate::lastD() const
  {
    return statsD_.lastValueMat();
  }

  DenseMatrix<double> ControlVariate::meanD() const
  {
    return statsD_.meanMat();
  }

  Vector<double> ControlVariate::lastGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.lastValueVec();
  }

  Vector<double> ControlVariate::meanGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.meanVec();
  }

  double ControlVariate::autocorrelationB2(int iOfSpan, int iOfFunction) const
  {
    return statsB2_.correlationAtSpan(iOfSpan, iOfFunction);
  }

  Vector<double> ControlVariate::lastA() const
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
    appendToB1(observable);
    appendToB2(observable, iOfStep);
    appendToD();
    appendToObservable(observable, iOfStep);
    appendToBetter(observable, iOfStep);
    //cout << "end ControlVariate::update" << endl;
  }
  


  void ControlVariate::display(std::ofstream& out, double time) const
  {
    out << time;
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanB()(iOfFunction) * lastA(iOfFunction);
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << lastA(iOfFunction);
    /*for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << lastB()(iOfFunction);*/
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanB()(iOfFunction);
    /*for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      for (int iOfFunction2 = 0; iOfFunction2 < nbOfFunctions(); iOfFunction2++)
    out << " " << lastD()(iOfFunction, iOfFunction2);*/
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      for (int iOfFunction2 = 0; iOfFunction2 < nbOfFunctions(); iOfFunction2++)
        out << " " << meanD()(iOfFunction, iOfFunction2);  //8-11

    out << " " << lastValue()
        << " " << mean()   //13
        << " " << variance()
        << " " << varOfVar()
        << " " << lastValueBetter()
        << " " << meanBetter()
        << " " << varOfVarBetter();
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << lastGeneratorOnBasis()(iOfFunction);
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanGeneratorOnBasis()(iOfFunction);

    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanB1()(iOfFunction);

    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << -correlationB2()(iOfFunction);  //22-23

    out << endl;
  }



 

  SinusControlVariate::SinusControlVariate(Input const& input, const string& outPath, CVBasis& cvBasis):
    ControlVariate(input, outPath, cvBasis, 1)
  {}

  double SinusControlVariate::basisFunction(vector<Particle*> const& configuration, int /*iOfFunction*/) const
  {
    double q = configuration[0]->position(0);
    return sin(2 * M_PI * q);
  }

  Vector<double> SinusControlVariate::gradientQ(vector<Particle*> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    double q = configuration[0]->position(0);
    Vector<double> grad(dimension_);
    grad(0) = 2 * M_PI * cos(2 * M_PI * q);
    return grad;
  }

  double SinusControlVariate::laplacianQ(vector<Particle*> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    double q = configuration[0]->position(0);
    return - pow (2 * M_PI, 2) * sin(2 * M_PI * q);
  }


  double SinusControlVariate::laplacianP(vector<Particle*> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return 0;
  }

  Vector<double> SinusControlVariate::gradientP(vector<Particle*> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return Vector<double>(dimension_);
  }


  //#####BasisControlVariate#####

  BasisControlVariate::BasisControlVariate(Input const& input, const string& outPath, CVBasis& cvBasis0, Galerkin* galerkin):
    ControlVariate(input, outPath, cvBasis0, 1),
    coeffsVec_(input.nbOfFourier() * input.nbOfHermite(), 1)
  {
    if (galerkin)
    {
      cout << "CVcoeffs from Galerkin !" << endl;
      coeffsVec_ = galerkin->CVcoeffs();
    }
    else
    {
      cout << "CVcoeffs from file !" << endl;
      std::string coeffsPath = input.outputFolderName() + input.controlVariateCoeffsPath();
      ifstream file(coeffsPath);
      coeffsVec_ = SparseMatrix<double>(coeffsPath, getNbOfLines(file));
    }
    //cout << "coeffsVec_ : " << endl;
    //cout << coeffsVec_ << endl;
  }

  double BasisControlVariate::basisFunction(vector<Particle*> const& configuration, int iOfFunction) const
  {
    assert(iOfFunction == 0);
    double result = 0;

    for (SMat::iterator it(coeffsVec_, 0); it; ++it)
    {
      int iOfCoeff = it.row();
      double valOfCoeff = it.value();
      result += valOfCoeff * cvBasis_->basis_->value(configuration, iOfCoeff);
    }

    return result;
  }

  Vector<double> BasisControlVariate::gradientQ(vector<Particle*> const& configuration, int iOfParticle, int iOfFunction) const
  {
    //cout << "BasisControlVariate::gradientQ" << endl;
    Vector<double> result(1, 0);
    assert(iOfFunction == 0);
    //cout << "gradientQ" << endl;
    for (SMat::iterator it(coeffsVec_, 0); it; ++it)
    {
      int iOfCoeff = it.row();
      double valOfCoeff = it.value();
      //cout << "+ " << valOfCoeff << " X ";
      result += valOfCoeff * cvBasis_->basis_->gradientQ(configuration, iOfParticle, iOfCoeff);
    }
    //cout << endl;
    //cout << "end BasisControlVariate::gradientQ" << endl;
    //cout << "gradientQ = " << result << endl;
    return result;
  }

  double BasisControlVariate::laplacianQ(vector<Particle*> const& configuration, int iOfParticle, int iOfFunction) const
  {
    double result = 0;
    assert(iOfFunction == 0);

    for (SMat::iterator it(coeffsVec_, 0); it; ++it)
    {
      int iOfCoeff = it.row();
      double valOfCoeff = it.value();
      result += valOfCoeff * cvBasis_->basis_->laplacianQ(configuration, iOfParticle, iOfCoeff);
    }
    //cout << "laplacianQ = " << result << endl;
    return result;
  }

  Vector<double> BasisControlVariate::gradientP(vector<Particle*> const& configuration, int iOfParticle, int iOfFunction) const
  {
    Vector<double> result(1, 0);
    assert(iOfFunction == 0);
    //cout << "gradientP" << endl;
    for (SMat::iterator it(coeffsVec_, 0); it; ++it)
    {
      int iOfCoeff = it.row();
      double valOfCoeff = it.value();
      //cout << "+ " << valOfCoeff << " X ";
      result += valOfCoeff * cvBasis_->basis_->gradientP(configuration, iOfParticle, iOfCoeff);
    }
    //cout << endl;
    //cout << "gradientP = " << result << endl;
    return result;
  }

  double BasisControlVariate::laplacianP(vector<Particle*> const& configuration, int iOfParticle, int iOfFunction) const
  {
    double result = 0;
    assert(iOfFunction == 0);

    for (SMat::iterator it(coeffsVec_, 0); it; ++it)
    {
      int iOfCoeff = it.row();
      double valOfCoeff = it.value();
      result += valOfCoeff * cvBasis_->basis_->laplacianP(configuration, iOfParticle, iOfCoeff);
    }
    //cout << "laplacianP = " << result << endl;
    return result;
  }



  ExpFourierHermiteControlVariate::ExpFourierHermiteControlVariate(Input const& input, const string& outPath, CVBasis& cvBasis0, Galerkin* galerkin):
    BasisControlVariate(input, outPath, cvBasis0, galerkin)
  {
    //basis_ = new ExpFourierHermiteBasis(input, potential);
    nbQ_ = 100;  // rafinement du maillage de l'output map
    nbP_ = 100;  // rafinement du maillage de l'output map
    deltaQ_ = 2 * M_PI / nbQ_;
    pMax_ = 4;
    deltaP_ = 2 * pMax_ / nbP_;
  }

  int ExpFourierHermiteControlVariate::nbOfFourier() const
  {return cvBasis_->basis_->nbOfElts(0);}

  int ExpFourierHermiteControlVariate::nbOfHermite() const
  {return cvBasis_->basis_->nbOfElts(1);}

  void ExpFourierHermiteControlVariate::displayMap(ofstream& out) const
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
        conf[0]->position() = Vector<double>(1, q);
        conf[0]->momentum() = Vector<double>(1, p);
        out << basisFunction(conf, 0) << " ";
      }
      out << endl;
    }
    cout << "end displayMap" << endl;
  }

  void ExpFourierHermiteControlVariate::displayGradQMap(ofstream& out) const
  {
    cout << "displayGradQMap" << endl;
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
        conf[0]->position() = Vector<double>(1, q);
        conf[0]->momentum() = Vector<double>(1, p);
        out << gradientQ(conf, 0) << " ";
      }
      out << endl;
    }
    cout << "end displayMap" << endl;
  }

  void ExpFourierHermiteControlVariate::displayGradPMap(ofstream& out) const
  {
    cout << "displayGradPMap" << endl;
    vector<Particle*> conf(1, new Particle(0, 0, 0));
    //nf[0] = Particle(0,0,0);
    out << nbP_ + 1 << " ";
    for (double p = -pMax_; p < pMax_; p += deltaP_)
      out << p << " ";
    out << endl;
    for (double q = -M_PI; q < M_PI; q += deltaQ_)
    {
      out << q << " ";
      for (double p = -pMax_; p < pMax_; p += deltaP_)
      {
        conf[0]->position() = Vector<double>(1, q);
        conf[0]->momentum() = Vector<double>(1, p);
        out << gradientP(conf, 0) << " ";
      }
      out << endl;
    }
    cout << "end displayMap" << endl;
  }

}










