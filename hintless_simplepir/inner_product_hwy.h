// Modified by Mayank Rathee
// Copyright 2024 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef HINTLESS_PIR_HINTLESS_SIMPLEPIR_INNER_PRODUCT_HWY_H_
#define HINTLESS_PIR_HINTLESS_SIMPLEPIR_INNER_PRODUCT_HWY_H_

#include <stdint.h>

#include <string>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "lwe/types.h"

namespace hintless_pir {
namespace hintless_simplepir {
namespace internal {

using BlockType = absl::uint128;
using BlockVector = std::vector<BlockType>;

// Given a matrix represented by its columns in `matrix`, and a vector `vec`,
// returns the product = `matrix` * `vec` (mod Q), where Q is the LWE modulus.
// The matrix stores its elements in PlainInteger (uint8_t or uint16_t), packed
// in BlockType; so each column is represented as a vector of BlockType.
// This version is implemented using SIMD instructions via the highway library.
template <typename PlainInteger>
absl::StatusOr<std::vector<lwe::Integer>> InnerProduct(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec, bool secondary = false);

template <typename PlainInteger>
absl::StatusOr<std::vector<lwe::Integer64>> InnerProduct64(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer64> vec, bool secondary = false);

/*
template <typename PlainInteger>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductFast32by8(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec_upper, absl::Span<const uint16_t> vec_lower);

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductHwyLazy(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec);

// For 54-bit prime modulus
template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductHwy_54(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec);
*/

// Matrix-vector product implemented without using highway SIMD intrinsics.
template <typename PlainInteger>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductNoHwyLazy(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec, bool secondary = false);

// For 54-bit prime modulus
template <typename PlainInteger>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductNoHwy_54(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec);

template <typename PlainInteger>
absl::StatusOr<std::vector<lwe::Integer64>> InnerProductNoHwy_54_64(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer64> vec);

template <typename PlainInteger>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductNoHwy_Native(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec);

/*
template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductNoHwyLazy32by8(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec_upper, absl::Span<const uint16_t> vec_lower);

template <typename PlainInteger>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductNoHwy_Solinas32Mod(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec);
*/

}  // namespace internal
}  // namespace hintless_simplepir
}  // namespace hintless_pir

#endif  // HINTLESS_PIR_HINTLESS_SIMPLEPIR_INNER_PRODUCT_HWY_H_
