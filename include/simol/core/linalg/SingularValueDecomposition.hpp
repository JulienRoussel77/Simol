#ifndef SIMOL_SINGULARVALUEDECOMPOSITION_HPP
#define SIMOL_SINGULARVALUEDECOMPOSITION_HPP

namespace simol
{
  template<typename ScalarType, template<class> class WrappedLibrary = eigen>
  class SingularValueDecomposition;

  template<typename ScalarType>
  class SingularValueDecomposition<ScalarType, eigen>
  {
    public:
      SingularValueDecomposition(DenseMatrix const & matrix);
    private:
      typename eigen<ScalarType>::SVDType wrapped_;
  };
}


#endif

