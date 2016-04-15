#include "simol/core/linalg/SparseTensor.hpp"

namespace simol
{
    std::size_t
    getInd(std::size_t const M_disc,
           std::size_t const rowIndex,
           std::size_t const columnIndex)
    {
        assert(rowIndex < M_disc && columnIndex < M_disc);
        return rowIndex * M_disc + columnIndex;
    }
}
