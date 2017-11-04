#include "simol/statphys/controlVariate/CVBasis.hpp"
#include <simol/core/io/CommandLine.hpp>

namespace simol
{
  shared_ptr<CVBasis> createCVBasis(Input const& input)
  {
    if (input.basisElts() == "None")
      return nullptr;
    return shared_ptr<CVBasis>(new CVBasis(input));
    
//     if (input.basisElts() == "None")
//       return nullptr;
//     else if (input.systemName() == "Colloid")
//       return make_shared<ColloidCVBasis>(input, dynaPara);
//     else if (input.systemName() == "Isolated")
//       return make_shared<IsolatedCVBasis>(input, dynaPara);
//     else 
//       throw runtime_error("createCVBasis does not know what CVBasis should be built !");
  }
  
  CVBasis::CVBasis(Input const& input):
    potential_(createPotential(input, input.galerkinPotentialName())),
    tensorBasis_(createTensorBasis(input, potential_)),
    totalNbOfElts_(tensorBasis_?tensorBasis_->totalNbOfElts():0),
    cvCoeffs_(nullptr),
    basisValues_(totalNbOfElts_),
    generatorOnBasisValues_(totalNbOfElts_)
  {
    generator_ = createCVOperator(input);
    cout << "Galerkin Potential : " << potential_->classname() << endl;
    if (input.controlVariateCoeffsPath() != "None")
    {
      if (input.basisElts() == "QuadraticHermite")
        cvCoeffs_ = make_shared<DVec>(DVec::Constant(1, 1, 1));
      else{
        cout << "CVcoeffs from file !" << endl;
        std::string coeffsPath = input.outputFolderName() + input.controlVariateCoeffsPath();
        vector<int> dimensions;
        cvCoeffs_ = make_shared<DVec>(scanTensor(coeffsPath, dimensions));
      }
    }
  }

  
  CVBasis::~CVBasis()
  {
    //if (basis_) delete basis_;
  }

  const int& CVBasis::totalNbOfElts() const
  {return totalNbOfElts_;}
  

  
  ///
  ///Evaluates each function of the basis at "configuration" and updates the attribute
  void CVBasis::computeValueBasis(System const& syst)
  {  
    DVec variables = generator_->basisVariables(syst);
    for (int iOfElt = 0; iOfElt < totalNbOfElts(); iOfElt++)
      basisValues_(iOfElt) = tensorBasis_->value(variables, iOfElt);
  }
  
  ///
  ///Evaluates the image of each function by the generator of the basis at "configuration" and updates the attribute
  void CVBasis::computeValueGeneratorOnBasis(System const& syst)
  {
    DVec variables = generator_->basisVariables(syst);
    generatorOnBasisValues_ = generator_->value(tensorBasis_, syst);
  }
  
  DMat CVBasis::gradientQ(System const& syst, int iOfFunction)
  {
    return tensorBasis_->gradientQ(generator_->basisVariables(syst), iOfFunction);    
  }
  
  DMat CVBasis::gradientP(System const& syst, int iOfFunction)
  {
    return tensorBasis_->gradientP(generator_->basisVariables(syst), iOfFunction);    
  }
  
  double CVBasis::laplacianQ(System const& syst, int iOfFunction)
  {
    return tensorBasis_->laplacianQ(generator_->basisVariables(syst), iOfFunction);
  }
  
  double CVBasis::laplacianP(System const& syst, int iOfFunction)
  {
    return tensorBasis_->laplacianP(generator_->basisVariables(syst), iOfFunction);    
  }
  
  
}   
  