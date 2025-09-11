#include "hintless_simplepir/matrix_utils.h"

#include <cstdint>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "lwe/types.h"
#include "lib/pir/mat.h"

namespace hintless_pir
{
    namespace hintless_simplepir
    {
        template <typename T>
        absl::StatusOr<Matrix> MatrixVectorTransposeMul(
            const Matrix &matrix,
            const Matrix &vec)
        {

            // Verify input is a column vector
            if (vec.cols != 1)
            {
                return absl::InvalidArgumentError("Second argument must be a column vector");
            }

            // Create result matrix (matrix.rows × vec.rows)
            Matrix result(matrix.rows, vec.rows);

            // Choose upcast type based on input type
            using UpcastType = typename std::conditional<
                std::is_same<T, uint64_t>::value,
                absl::uint128_t,
                uint64_t>::type;

            // Compute r^T * Matrix (transpose vector times matrix)
            for (size_t i = 0; i < matrix.rows; i++)
            {
                for (size_t j = 0; j < vec.rows; j++)
                {
                    // Upcast both operands before multiplication
                    UpcastType matrix_val = static_cast<UpcastType>(matrix.data[matrix.cols * i]);
                    UpcastType vec_val = static_cast<UpcastType>(vec.data[j]);

                    // Multiply upcasted values and downcast result
                    result.data[i * vec.rows + j] = static_cast<Elem>(fast_mod(matrix_val * vec_val));
                }
            }

            return result;
        }

        // Explicit template instantiations
        template absl::StatusOr<Matrix> MatrixVectorTransposeMul<uint32_t>(
            const Matrix &matrix,
            const Matrix &vec);

        template absl::StatusOr<Matrix> MatrixVectorTransposeMul<uint64_t>(
            const Matrix &matrix,
            const Matrix &vec);

    } // namespace hintless_simplepir
} // namespace hintless_pir