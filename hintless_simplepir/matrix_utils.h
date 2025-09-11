#ifndef HINTLESS_SIMPLEPIR_SRC_LIB_PIR_MATRIX_UTILS_H_
#define HINTLESS_SIMPLEPIR_SRC_LIB_PIR_MATRIX_UTILS_H_

#include <cstdint>
#include "absl/status/statusor.h"
#include "lib/pir/mat.h"

namespace hintless_pir
{
    namespace hintless_simplepir
    {

        // Computes r^T * matrix, where r^T is the transpose of vector r.
        // Uses 128-bit integers for intermediate values when input type is uint64_t,
        // and 64-bit integers when input type is uint32_t.
        //
        // Template parameter T must be either uint32_t or uint64_t.
        //
        // Returns:
        //   A StatusOr containing the result matrix on success, or an error Status if:
        //   - Input vector is not a column vector (vec.cols != 1)
        template <typename T>
        absl::StatusOr<Matrix> MatrixVectorTransposeMul(
            const Matrix &matrix,
            const Matrix &vec);

    } // namespace hintless_simplepir
} // namespace hintless_pir

#endif // HINTLESS_SIMPLEPIR_MATRIX_VECTOR_TRANSPOSE_H_