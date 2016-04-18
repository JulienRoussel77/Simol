#ifndef SIMOL_ARPACK_HPP
#define SIMOL_ARPACK_HPP

#include "dsaupd.hpp"
#include "dseupd.hpp"

#include "EigenSolver.hpp"

namespace simol
{
  class Arpack
    : public EigenSolver
  {
    public:
      Arpack();
    private:
      int ncv;

  };
}


#endif
