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

#include "hintless_simplepir/inner_product_hwy.h"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "absl/base/optimization.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "hwy/detect_targets.h"
#include "lwe/types.h"

// Highway implementations.
// clang-format off
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "hintless_simplepir/inner_product_hwy.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
// clang-format on

// Must come after foreach_target.h to avoid redefinition errors.
#include "hwy/aligned_allocator.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace hintless_pir::hintless_simplepir::internal {
namespace HWY_NAMESPACE {

#if HWY_TARGET == HWY_SCALAR

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductHwy(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec, bool secondary) {
      if (SOL_NTT_PRIME_BITS == 54) {
        return InnerProductNoHwy_54<DBIntegerType>(matrix, vec);  
      }

      if (SOL_NTT_PRIME_BITS == 0) {
        return InnerProductNoHwy_Native<DBIntegerType>(matrix, vec);  
      }

      uint64_t max = -1;
      //assert(max > (vec.size() * SOL_NTT_PRIME));
      if((1ULL<<(8*sizeof(DBIntegerType))) >= max/(vec.size() * (SOL_NTT_PRIME))) {
        std::cout << "\n ERROR: Inner product correctness compromised from lazy reduction.\n" << std::endl;
      }
      
      return InnerProductNoHwyLazy<DBIntegerType>(matrix, vec, secondary);
}

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer64>> InnerProductHwy64(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer64> vec, bool secondary) {
      if (VSPIR_LHE_PRIME_BITS != 54) {
        std::cout << "\n ERROR: Inner product 64-bit only supported for 54-bit prime modulus.\n" << std::endl;
        assert(false);  
      }

      return InnerProductNoHwy_54_64<DBIntegerType>(matrix, vec);
}

#else

namespace hn = hwy::HWY_NAMESPACE;

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductHwy_Native(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec) {
  if (matrix.size() != vec.size()) {
    return absl::InvalidArgumentError(
        "`matrix` and `vec` must have matching dimensions.");
  }

  // Vector type used throughout this function: Largest byte vector
  // available.
  const hn::ScalableTag<lwe::Integer> d32;
  const hn::Rebind<DBIntegerType, hn::ScalableTag<lwe::Integer>> d_plain;
  const int N = hn::Lanes(d32);

  // Do not run the highway version if
  // - the number of bytes in a hwy vector is less than 16, or
  // - the number of bytes in a hwy vector is not a multiple of 16.
  if (ABSL_PREDICT_FALSE(N < 4 || N % 4 != 0)) {
    return InnerProductNoHwy_Native<DBIntegerType>(matrix, vec);
  }

  // Assume all columns have the same size.
  int num_blocks = matrix[0].size();
  int num_values_per_block = sizeof(BlockType) / sizeof(DBIntegerType);
  int num_rows = num_blocks * num_values_per_block;

  // Allocate aligned buffers to hold the intermediate values of inner
  // products.
  hwy::AlignedFreeUniquePtr<lwe::Integer[]> aligned_results =
      hwy::AllocateAligned<lwe::Integer>(num_rows);
  std::fill_n(aligned_results.get(), num_rows, 0);

  for (int j = 0; j < vec.size(); ++j) {
    int row_idx = 0;
    const DBIntegerType* value_ptr = 0;
    // First, run 4x SIMD multiplication in each iteration.
    for (; row_idx + N * 4 <= num_rows; row_idx += N * 4) {
      if (row_idx % num_values_per_block == 0) {
        int block_idx = row_idx / num_values_per_block;
        value_ptr =
            reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      } else {
        value_ptr += N;
      }
      lwe::Integer* result_ptr = &aligned_results[row_idx];
      auto add32_0 = hn::Load(d32, result_ptr);
      auto add32_1 = hn::Load(d32, result_ptr + N);
      auto add32_2 = hn::Load(d32, result_ptr + 2 * N);
      auto add32_3 = hn::Load(d32, result_ptr + 3 * N);

      auto left0 = hn::LoadU(d_plain, value_ptr);
      auto left1 = hn::LoadU(d_plain, value_ptr + N);
      auto left2 = hn::LoadU(d_plain, value_ptr + 2 * N);
      auto left3 = hn::LoadU(d_plain, value_ptr + 3 * N);

      auto left32_0 = hn::PromoteTo(d32, left0);
      auto left32_1 = hn::PromoteTo(d32, left1);
      auto left32_2 = hn::PromoteTo(d32, left2);
      auto left32_3 = hn::PromoteTo(d32, left3);

      auto right32 = hn::Set(d32, vec[j]);

      auto mul32_0 = hn::MulAdd(left32_0, right32, add32_0);
      auto mul32_1 = hn::MulAdd(left32_1, right32, add32_1);
      auto mul32_2 = hn::MulAdd(left32_2, right32, add32_2);
      auto mul32_3 = hn::MulAdd(left32_3, right32, add32_3);

      hn::Store(mul32_0, d32, result_ptr);
      hn::Store(mul32_1, d32, result_ptr + N);
      hn::Store(mul32_2, d32, result_ptr + 2 * N);
      hn::Store(mul32_3, d32, result_ptr + 3 * N);
    }

    // Next, run 1x per iteration.
    for (; row_idx + N <= num_rows; row_idx += N) {
      if (row_idx % num_values_per_block == 0) {
        int block_idx = row_idx / num_values_per_block;
        value_ptr =
            reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      } else {
        value_ptr += N;
      }
      lwe::Integer* result_ptr = &aligned_results[row_idx];
      auto add32 = hn::Load(d32, result_ptr);
      auto left = hn::LoadU(d_plain, value_ptr);
      auto left32 = hn::PromoteTo(d32, left);
      auto right32 = hn::Set(d32, vec[j]);
      auto mul32 = hn::MulAdd(left32, right32, add32);
      hn::Store(mul32, d32, result_ptr);
    }

    // Handle the remaining rows that didn't take a full lane.
    if (row_idx < num_rows) {
      int block_idx = row_idx / num_values_per_block;
      const DBIntegerType* block_as_values =
          reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      for (; row_idx < num_rows; ++row_idx) {
        if (row_idx % num_values_per_block == 0) {
          // update the block pointer, which should be rate
          block_idx = row_idx / num_values_per_block;
          block_as_values =
              reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
        }
        int block_pos = row_idx % num_values_per_block;
        aligned_results.get()[row_idx] +=
            static_cast<lwe::Integer>(block_as_values[block_pos]) * vec[j];
      }
    }
  }

  return std::vector<lwe::Integer>(aligned_results.get(),
                                   aligned_results.get() + num_rows);
}

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductHwyLazy(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec, bool secondary = false) {

  std::function<uint64_t(uint64_t)> mod_func_max_64 = fast_mod_max_64;
  if (secondary) {
    mod_func_max_64 = fast_mod_secondary_max_64;
  }

  if (matrix.size() != vec.size()) {
    return absl::InvalidArgumentError(
        "`matrix` and `vec` must have matching dimensions.");
  }

  // For lazy reduction correctness
  uint64_t max = -1;
  assert(max > (vec.size() * SOL_NTT_PRIME));
  assert((1ULL<<(8*sizeof(DBIntegerType))) < max/(vec.size() * (SOL_NTT_PRIME)));

  // Vector type used throughout this function: Largest byte vector
  // available.
  const hn::ScalableTag<Elem64> d64;
  const hn::Rebind<DBIntegerType, hn::ScalableTag<Elem64>> d_plain;
  const int N = hn::Lanes(d64);

  // Do not run the highway version if
  // - the number of bytes in a hwy vector is less than 16, or
  // - the number of bytes in a hwy vector is not a multiple of 16.
  if (ABSL_PREDICT_FALSE(N < 4 || N % 4 != 0)) {
    return InnerProductNoHwyLazy<DBIntegerType>(matrix, vec);
  }

  // Assume all columns have the same size.
  int num_blocks = matrix[0].size();
  int num_values_per_block = sizeof(BlockType) / sizeof(DBIntegerType);
  int num_rows = num_blocks * num_values_per_block;

  std::vector<lwe::Integer> result(num_rows, 0);

  // Allocate aligned buffers to hold the intermediate values of inner
  // products.
  hwy::AlignedFreeUniquePtr<Elem64[]> aligned_results =
      hwy::AllocateAligned<Elem64>(num_rows);
  std::fill_n(aligned_results.get(), num_rows, 0);

  for (int j = 0; j < vec.size(); ++j) {
    int row_idx = 0;
    const DBIntegerType* value_ptr = 0;
    // First, run 4x SIMD multiplication in each iteration.
    for (; row_idx + N * 4 <= num_rows; row_idx += N * 4) {
      if (row_idx % num_values_per_block == 0) {
        int block_idx = row_idx / num_values_per_block;
        value_ptr =
            reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      } else {
        value_ptr += N;
      }
      Elem64* result_ptr = &aligned_results[row_idx];
      auto add64_0 = hn::Load(d64, result_ptr);
      auto add64_1 = hn::Load(d64, result_ptr + N);
      auto add64_2 = hn::Load(d64, result_ptr + 2 * N);
      auto add64_3 = hn::Load(d64, result_ptr + 3 * N);

      auto left0 = hn::LoadU(d_plain, value_ptr);
      auto left1 = hn::LoadU(d_plain, value_ptr + N);
      auto left2 = hn::LoadU(d_plain, value_ptr + 2 * N);
      auto left3 = hn::LoadU(d_plain, value_ptr + 3 * N);

      auto left64_0 = hn::PromoteTo(d64, left0);
      auto left64_1 = hn::PromoteTo(d64, left1);
      auto left64_2 = hn::PromoteTo(d64, left2);
      auto left64_3 = hn::PromoteTo(d64, left3);

      auto right64 = hn::Set(d64, vec[j]);

      auto mul64_0 = hn::MulAdd(left64_0, right64, add64_0);
      auto mul64_1 = hn::MulAdd(left64_1, right64, add64_1);
      auto mul64_2 = hn::MulAdd(left64_2, right64, add64_2);
      auto mul64_3 = hn::MulAdd(left64_3, right64, add64_3);

      hn::Store(mul64_0, d64, result_ptr);
      hn::Store(mul64_1, d64, result_ptr + N);
      hn::Store(mul64_2, d64, result_ptr + 2 * N);
      hn::Store(mul64_3, d64, result_ptr + 3 * N);
    }

    // Next, run 1x per iteration.
    for (; row_idx + N <= num_rows; row_idx += N) {
      if (row_idx % num_values_per_block == 0) {
        int block_idx = row_idx / num_values_per_block;
        value_ptr =
            reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      } else {
        value_ptr += N;
      }
      Elem64* result_ptr = &aligned_results[row_idx];
      auto add64 = hn::Load(d64, result_ptr);
      auto left = hn::LoadU(d_plain, value_ptr);
      auto left64 = hn::PromoteTo(d64, left);
      auto right64 = hn::Set(d64, vec[j]);
      auto mul64 = hn::MulAdd(left64, right64, add64);
      hn::Store(mul64, d64, result_ptr);
    }

    // Handle the remaining rows that didn't take a full lane.
    if (row_idx < num_rows) {
      int block_idx = row_idx / num_values_per_block;
      const DBIntegerType* block_as_values =
          reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      for (; row_idx < num_rows; ++row_idx) {
        if (row_idx % num_values_per_block == 0) {
          // update the block pointer, which should be rate
          block_idx = row_idx / num_values_per_block;
          block_as_values =
              reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
        }
        int block_pos = row_idx % num_values_per_block;
        aligned_results.get()[row_idx] +=
            static_cast<Elem64>(block_as_values[block_pos]) * vec[j];
      }
    }
  }

  auto result_unreduced = std::vector<Elem64>(aligned_results.get(),
                                   aligned_results.get() + num_rows);

  for (int i = 0; i < num_rows; i++){
    result[i] = static_cast<lwe::Integer>(mod_func_max_64(result_unreduced[i]));
  }
  return result;
}

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductHwy_54(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec) {

  if (matrix.size() != vec.size()) {
    return absl::InvalidArgumentError(
        "`matrix` and `vec` must have matching dimensions.");
  }

  // Vector type used throughout this function: Largest byte vector
  // available.
  const hn::ScalableTag<Elem64> d64;
  const hn::Rebind<DBIntegerType, hn::ScalableTag<Elem64>> d_plain;
  const int N = hn::Lanes(d64);

  // Do not run the highway version if
  // - the number of bytes in a hwy vector is less than 16, or
  // - the number of bytes in a hwy vector is not a multiple of 16.
  if (ABSL_PREDICT_FALSE(N < 4 || N % 4 != 0)) {
    return InnerProductNoHwy_54<DBIntegerType>(matrix, vec);
  }

  // Assume all columns have the same size.
  int num_blocks = matrix[0].size();
  int num_values_per_block = sizeof(BlockType) / sizeof(DBIntegerType);
  int num_rows = num_blocks * num_values_per_block;

  std::vector<lwe::Integer> result(num_rows, 0);

  // Allocate aligned buffers to hold the intermediate values of inner
  // products.
  hwy::AlignedFreeUniquePtr<Elem64[]> aligned_results =
      hwy::AllocateAligned<Elem64>(num_rows);
  std::fill_n(aligned_results.get(), num_rows, 0);
  Elem64 tmp_var = 0;
  Elem64 mask_54 = (1ULL << 54) - 1;

  for (int j = 0; j < vec.size(); ++j) {
    int row_idx = 0;
    const DBIntegerType* value_ptr = 0;
    // First, run 4x SIMD multiplication in each iteration.
    for (; row_idx + N * 4 <= num_rows; row_idx += N * 4) {
      if (row_idx % num_values_per_block == 0) {
        int block_idx = row_idx / num_values_per_block;
        value_ptr =
            reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      } else {
        value_ptr += N;
      }
      Elem64* result_ptr = &aligned_results[row_idx];
      auto add64_0 = hn::Load(d64, result_ptr);
      auto add64_1 = hn::Load(d64, result_ptr + N);
      auto add64_2 = hn::Load(d64, result_ptr + 2 * N);
      auto add64_3 = hn::Load(d64, result_ptr + 3 * N);

      auto left0 = hn::LoadU(d_plain, value_ptr);
      auto left1 = hn::LoadU(d_plain, value_ptr + N);
      auto left2 = hn::LoadU(d_plain, value_ptr + 2 * N);
      auto left3 = hn::LoadU(d_plain, value_ptr + 3 * N);

      auto left64_0 = hn::PromoteTo(d64, left0);
      auto left64_1 = hn::PromoteTo(d64, left1);
      auto left64_2 = hn::PromoteTo(d64, left2);
      auto left64_3 = hn::PromoteTo(d64, left3);

      auto right64 = hn::Set(d64, vec[j]);

      auto mul64_0 = hn::MulAdd(left64_0, right64, add64_0);
      auto mul64_1 = hn::MulAdd(left64_1, right64, add64_1);
      auto mul64_2 = hn::MulAdd(left64_2, right64, add64_2);
      auto mul64_3 = hn::MulAdd(left64_3, right64, add64_3);

      // Inlined modulo reduction without any comparisons
      auto mask_54_bcast = hn::Set(d64, mask_54);

      auto lower_0 = mul64_0 & mask_54_bcast;
      auto lower_1 = mul64_1 & mask_54_bcast;
      auto lower_2 = mul64_2 & mask_54_bcast;
      auto lower_3 = mul64_3 & mask_54_bcast;

      auto upper_0 = hn::ShiftRightSame(mul64_0, 54);
      auto upper_1 = hn::ShiftRightSame(mul64_1, 54);
      auto upper_2 = hn::ShiftRightSame(mul64_2, 54);
      auto upper_3 = hn::ShiftRightSame(mul64_3, 54);

      auto upper_24_0 = hn::ShiftLeftSame(upper_0, 24);
      auto upper_24_1 = hn::ShiftLeftSame(upper_1, 24);
      auto upper_24_2 = hn::ShiftLeftSame(upper_2, 24);
      auto upper_24_3 = hn::ShiftLeftSame(upper_3, 24);

      auto reduced_0 = lower_0 - upper_0 + upper_24_0;
      auto reduced_1 = lower_1 - upper_1 + upper_24_1;
      auto reduced_2 = lower_2 - upper_2 + upper_24_2;
      auto reduced_3 = lower_3 - upper_3 + upper_24_3; 

      hn::Store(reduced_0, d64, result_ptr);
      hn::Store(reduced_1, d64, result_ptr + N);
      hn::Store(reduced_2, d64, result_ptr + 2 * N);
      hn::Store(reduced_3, d64, result_ptr + 3 * N);
    }

    // Next, run 1x per iteration.
    for (; row_idx + N <= num_rows; row_idx += N) {
      if (row_idx % num_values_per_block == 0) {
        int block_idx = row_idx / num_values_per_block;
        value_ptr =
            reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      } else {
        value_ptr += N;
      }
      Elem64* result_ptr = &aligned_results[row_idx];
      auto add64 = hn::Load(d64, result_ptr);
      auto left = hn::LoadU(d_plain, value_ptr);
      auto left64 = hn::PromoteTo(d64, left);
      auto right64 = hn::Set(d64, vec[j]);
      auto mul64 = hn::MulAdd(left64, right64, add64);
      // Inlined modulo reduction without any comparisons
      auto mask_54_bcast = hn::Set(d64, mask_54);
      auto lower = mul64 & mask_54_bcast;
      auto upper = hn::ShiftRightSame(mul64, 54);
      auto upper_24 = hn::ShiftLeftSame(upper, 24);
      auto reduced = lower - upper + upper_24;

      hn::Store(reduced, d64, result_ptr);
    }

    // Handle the remaining rows that didn't take a full lane.
    if (row_idx < num_rows) {
      int block_idx = row_idx / num_values_per_block;
      const DBIntegerType* block_as_values =
          reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      for (; row_idx < num_rows; ++row_idx) {
        if (row_idx % num_values_per_block == 0) {
          // update the block pointer, which should be rate
          block_idx = row_idx / num_values_per_block;
          block_as_values =
              reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
        }
        int block_pos = row_idx % num_values_per_block;
        //aligned_results.get()[row_idx] +=
        //    static_cast<Elem64>(block_as_values[block_pos]) * vec[j];
        tmp_var = aligned_results.get()[row_idx] + (static_cast<Elem64>(block_as_values[block_pos]) * vec[j]);
        // Inlined modulo reduction without any comparisons
        Elem64 lower = tmp_var & mask_54;
        Elem64 upper = tmp_var >> 54;
        aligned_results.get()[row_idx] = lower - upper + (upper << 24);

      }
    }
  }

  auto result_unreduced = std::vector<Elem64>(aligned_results.get(),
                                   aligned_results.get() + num_rows);

  for (int i = 0; i < num_rows; i++){
    result[i] = static_cast<lwe::Integer>(fast_mod_max_64(result_unreduced[i]));
  }
  return result;
}

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer64>> InnerProductHwy_54_64(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer64> vec) {

  if (matrix.size() != vec.size()) {
    return absl::InvalidArgumentError(
        "`matrix` and `vec` must have matching dimensions.");
  }

  // Vector type used throughout this function: Largest byte vector
  // available.
  const hn::ScalableTag<Elem64> d64;
  const hn::Rebind<DBIntegerType, hn::ScalableTag<Elem64>> d_plain;
  const int N = hn::Lanes(d64);

  // Do not run the highway version if
  // - the number of bytes in a hwy vector is less than 16, or
  // - the number of bytes in a hwy vector is not a multiple of 16.
  if (ABSL_PREDICT_FALSE(N < 4 || N % 4 != 0)) {
    return InnerProductNoHwy_54_64<DBIntegerType>(matrix, vec);
  }

  // Assume all columns have the same size.
  int num_blocks = matrix[0].size();
  int num_values_per_block = sizeof(BlockType) / sizeof(DBIntegerType);
  int num_rows = num_blocks * num_values_per_block;

  std::vector<lwe::Integer64> result(num_rows, 0);

  // Allocate aligned buffers to hold the intermediate values of inner
  // products.
  hwy::AlignedFreeUniquePtr<Elem64[]> aligned_results =
      hwy::AllocateAligned<Elem64>(num_rows);
  std::fill_n(aligned_results.get(), num_rows, 0);
  Elem64 tmp_var = 0;
  Elem64 mask_54 = (1ULL << 54) - 1;

  for (int j = 0; j < vec.size(); ++j) {
    int row_idx = 0;
    const DBIntegerType* value_ptr = 0;
    // First, run 4x SIMD multiplication in each iteration.
    for (; row_idx + N * 4 <= num_rows; row_idx += N * 4) {
      if (row_idx % num_values_per_block == 0) {
        int block_idx = row_idx / num_values_per_block;
        value_ptr =
            reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      } else {
        value_ptr += N;
      }
      Elem64* result_ptr = &aligned_results[row_idx];
      auto add64_0 = hn::Load(d64, result_ptr);
      auto add64_1 = hn::Load(d64, result_ptr + N);
      auto add64_2 = hn::Load(d64, result_ptr + 2 * N);
      auto add64_3 = hn::Load(d64, result_ptr + 3 * N);

      auto left0 = hn::LoadU(d_plain, value_ptr);
      auto left1 = hn::LoadU(d_plain, value_ptr + N);
      auto left2 = hn::LoadU(d_plain, value_ptr + 2 * N);
      auto left3 = hn::LoadU(d_plain, value_ptr + 3 * N);

      auto left64_0 = hn::PromoteTo(d64, left0);
      auto left64_1 = hn::PromoteTo(d64, left1);
      auto left64_2 = hn::PromoteTo(d64, left2);
      auto left64_3 = hn::PromoteTo(d64, left3);

      auto right64 = hn::Set(d64, vec[j]);

      auto mul64_0 = hn::MulAdd(left64_0, right64, add64_0);
      auto mul64_1 = hn::MulAdd(left64_1, right64, add64_1);
      auto mul64_2 = hn::MulAdd(left64_2, right64, add64_2);
      auto mul64_3 = hn::MulAdd(left64_3, right64, add64_3);

      // Inlined modulo reduction without any comparisons
      auto mask_54_bcast = hn::Set(d64, mask_54);

      auto lower_0 = mul64_0 & mask_54_bcast;
      auto lower_1 = mul64_1 & mask_54_bcast;
      auto lower_2 = mul64_2 & mask_54_bcast;
      auto lower_3 = mul64_3 & mask_54_bcast;

      auto upper_0 = hn::ShiftRightSame(mul64_0, 54);
      auto upper_1 = hn::ShiftRightSame(mul64_1, 54);
      auto upper_2 = hn::ShiftRightSame(mul64_2, 54);
      auto upper_3 = hn::ShiftRightSame(mul64_3, 54);

      auto upper_24_0 = hn::ShiftLeftSame(upper_0, 24);
      auto upper_24_1 = hn::ShiftLeftSame(upper_1, 24);
      auto upper_24_2 = hn::ShiftLeftSame(upper_2, 24);
      auto upper_24_3 = hn::ShiftLeftSame(upper_3, 24);

      auto reduced_0 = lower_0 - upper_0 + upper_24_0;
      auto reduced_1 = lower_1 - upper_1 + upper_24_1;
      auto reduced_2 = lower_2 - upper_2 + upper_24_2;
      auto reduced_3 = lower_3 - upper_3 + upper_24_3; 

      hn::Store(reduced_0, d64, result_ptr);
      hn::Store(reduced_1, d64, result_ptr + N);
      hn::Store(reduced_2, d64, result_ptr + 2 * N);
      hn::Store(reduced_3, d64, result_ptr + 3 * N);
    }

    // Next, run 1x per iteration.
    for (; row_idx + N <= num_rows; row_idx += N) {
      if (row_idx % num_values_per_block == 0) {
        int block_idx = row_idx / num_values_per_block;
        value_ptr =
            reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      } else {
        value_ptr += N;
      }
      Elem64* result_ptr = &aligned_results[row_idx];
      auto add64 = hn::Load(d64, result_ptr);
      auto left = hn::LoadU(d_plain, value_ptr);
      auto left64 = hn::PromoteTo(d64, left);
      auto right64 = hn::Set(d64, vec[j]);
      auto mul64 = hn::MulAdd(left64, right64, add64);
      // Inlined modulo reduction without any comparisons
      auto mask_54_bcast = hn::Set(d64, mask_54);
      auto lower = mul64 & mask_54_bcast;
      auto upper = hn::ShiftRightSame(mul64, 54);
      auto upper_24 = hn::ShiftLeftSame(upper, 24);
      auto reduced = lower - upper + upper_24;

      hn::Store(reduced, d64, result_ptr);
    }

    // Handle the remaining rows that didn't take a full lane.
    if (row_idx < num_rows) {
      int block_idx = row_idx / num_values_per_block;
      const DBIntegerType* block_as_values =
          reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
      for (; row_idx < num_rows; ++row_idx) {
        if (row_idx % num_values_per_block == 0) {
          // update the block pointer, which should be rate
          block_idx = row_idx / num_values_per_block;
          block_as_values =
              reinterpret_cast<const DBIntegerType*>(&matrix[j][block_idx]);
        }
        int block_pos = row_idx % num_values_per_block;
        //aligned_results.get()[row_idx] +=
        //    static_cast<Elem64>(block_as_values[block_pos]) * vec[j];
        tmp_var = aligned_results.get()[row_idx] + (static_cast<Elem64>(block_as_values[block_pos]) * vec[j]);
        // Inlined modulo reduction without any comparisons
        Elem64 lower = tmp_var & mask_54;
        Elem64 upper = tmp_var >> 54;
        aligned_results.get()[row_idx] = lower - upper + (upper << 24);

      }
    }
  }

  auto result_unreduced = std::vector<Elem64>(aligned_results.get(),
                                   aligned_results.get() + num_rows);

  for (int i = 0; i < num_rows; i++){
    result[i] = static_cast<lwe::Integer64>(fast_mod_vspir_max_64(result_unreduced[i]));
  }
  
  return result;
}

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductHwy(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec, bool secondary = false) {
      if (SOL_NTT_PRIME_BITS == 54) {
        return InnerProductHwy_54<DBIntegerType>(matrix, vec);  
      }

      if (SOL_NTT_PRIME_BITS == 0) {
        return InnerProductHwy_Native<DBIntegerType>(matrix, vec);  
      }

      uint64_t max = -1;
      //assert(max > (vec.size() * SOL_NTT_PRIME));
      if((1ULL<<(8*sizeof(DBIntegerType))) >= max/(vec.size() * (SOL_NTT_PRIME))) {
        std::cout << "\n ERROR: Inner product correctness compromised from lazy reduction.\n" << std::endl;
      }
      
      return InnerProductHwyLazy<DBIntegerType>(matrix, vec, secondary);
}

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer64>> InnerProductHwy64(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer64> vec, bool secondary = false) {
      if (VSPIR_LHE_PRIME_BITS != 54) {
        std::cout << "\n ERROR: Inner product 64-bit only supported for 54-bit prime modulus.\n" << std::endl;
        assert(false);
      }
      return InnerProductHwy_54_64<DBIntegerType>(matrix, vec);
}

#endif  // HWY_TARGET == HWY_SCALAR

}  // namespace HWY_NAMESPACE
}  // namespace hintless_pir::hintless_simplepir::internal
HWY_AFTER_NAMESPACE();

#if HWY_ONCE || HWY_IDE
namespace hintless_pir::hintless_simplepir::internal {

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductNoHwyLazy(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec, bool secondary) {

  std::function<uint64_t(uint64_t)> mod_func_max_64 = fast_mod_max_64;
  if (secondary) {
    mod_func_max_64 = fast_mod_secondary_max_64;
  }

  if (matrix.size() != vec.size()) {
    return absl::InvalidArgumentError(
        "`matrix` and `vec` must have matching dimensions.");
  }

  constexpr int num_values_per_block =
      sizeof(BlockType) / sizeof(DBIntegerType);

  // Assume all columns have the same size.
  int num_blocks = matrix[0].size();

  int num_rows = num_blocks * num_values_per_block;

  std::vector<Elem64> result_unreduced(num_rows, 0);
  std::vector<lwe::Integer> result(num_rows, 0);
  //Elem64 big_temp;
  for (int j = 0; j < vec.size(); ++j) {
    int i = 0;
    for (int block_idx = 0; block_idx < num_blocks; ++block_idx) {
      BlockType block = matrix[j][block_idx];
      const DBIntegerType* block_as_values =
          reinterpret_cast<const DBIntegerType*>(&block);
      // NOTE: This reinterpret cast is very unstable. It only really words for u8. For u16, it will not work, but we know that our database values are at most 8-bits, so we can use this for now.
      for (int block_pos = 0; block_pos < num_values_per_block && i < num_rows;
           ++block_pos, ++i) {
            result_unreduced[i] +=
              (static_cast<Elem64>(block_as_values[block_pos])) * vec[j];
      }
    }
  }

  for (int i = 0; i < num_rows; i++){
    result[i] = static_cast<lwe::Integer>(mod_func_max_64(result_unreduced[i]));
  }
  
  return result;
}

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductNoHwy_54(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec) {
  
  if (matrix.size() != vec.size()) {
    return absl::InvalidArgumentError(
        "`matrix` and `vec` must have matching dimensions.");
  }

  constexpr int num_values_per_block =
      sizeof(BlockType) / sizeof(DBIntegerType);

  // Assume all columns have the same size.
  int num_blocks = matrix[0].size();
  
  int num_rows = num_blocks * num_values_per_block;
  
  Elem64 mask_54 = (1ULL << 54) - 1;

  std::vector<Elem64> result_unreduced(num_rows, 0);
  std::vector<lwe::Integer> result(num_rows, 0);
  //Elem64 big_temp;
  for (int j = 0; j < vec.size(); ++j) {
    int i = 0;
    for (int block_idx = 0; block_idx < num_blocks; ++block_idx) {
      BlockType block = matrix[j][block_idx];
      const DBIntegerType* block_as_values =
          reinterpret_cast<const DBIntegerType*>(&block);
      // NOTE: This reinterpret cast is very unstable. It only really words for u8. For u16, it will not work, but we know that our database values are at most 8-bits, so we can use this for now.
      for (int block_pos = 0; block_pos < num_values_per_block && i < num_rows;
           ++block_pos, ++i) {
            result_unreduced[i] +=
              (static_cast<Elem64>(block_as_values[block_pos])) * vec[j];
            // Inlined modulo reduction without any comparisons
            Elem64 lower = result_unreduced[i] & mask_54;
            Elem64 upper = result_unreduced[i] >> 54;
            result_unreduced[i] = lower - upper + (upper << 24);
      }
    }
  }

  for (int i = 0; i < num_rows; i++){
    result[i] = static_cast<lwe::Integer>(fast_mod_max_64(result_unreduced[i]));
  }
  
  return result;
}

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer64>> InnerProductNoHwy_54_64(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer64> vec) {
  
  if (matrix.size() != vec.size()) {
    return absl::InvalidArgumentError(
        "`matrix` and `vec` must have matching dimensions.");
  }

  constexpr int num_values_per_block =
      sizeof(BlockType) / sizeof(DBIntegerType);

  // Assume all columns have the same size.
  int num_blocks = matrix[0].size();
  
  int num_rows = num_blocks * num_values_per_block;

  Elem64 mask_54 = (1ULL << 54) - 1;

  std::vector<Elem64> result_unreduced(num_rows, 0);
  std::vector<lwe::Integer64> result(num_rows, 0);
  //Elem64 big_temp;
  for (int j = 0; j < vec.size(); ++j) {
    int i = 0;
    for (int block_idx = 0; block_idx < num_blocks; ++block_idx) {
      BlockType block = matrix[j][block_idx];
      const DBIntegerType* block_as_values =
          reinterpret_cast<const DBIntegerType*>(&block);
      // NOTE: This reinterpret cast is very unstable. It only really words for u8. For u16, it will not work, but we know that our database values are at most 8-bits, so we can use this for now.
      for (int block_pos = 0; block_pos < num_values_per_block && i < num_rows;
           ++block_pos, ++i) {
            result_unreduced[i] +=
              (static_cast<Elem64>(block_as_values[block_pos])) * vec[j];
            // Inlined modulo reduction without any comparisons
            Elem64 lower = result_unreduced[i] & mask_54;
            Elem64 upper = result_unreduced[i] >> 54;
            result_unreduced[i] = lower - upper + (upper << 24);
      }
    }
  }

  for (int i = 0; i < num_rows; i++){
    result[i] = static_cast<lwe::Integer64>(fast_mod_vspir_max_64(result_unreduced[i]));
    //result[i] = static_cast<lwe::Integer>(result_unreduced[i]);
  }
  
  return result;
}

template <typename DBIntegerType>
absl::StatusOr<std::vector<lwe::Integer>> InnerProductNoHwy_Native(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec) {
  if (matrix.size() != vec.size()) {
    return absl::InvalidArgumentError(
        "`matrix` and `vec` must have matching dimensions.");
  }

  constexpr int num_values_per_block =
      sizeof(BlockType) / sizeof(DBIntegerType);

  // Assume all columns have the same size.
  int num_blocks = matrix[0].size();
  //int num_rows = num_blocks * sizeof(BlockType);
  int num_rows = num_blocks * num_values_per_block;

  std::vector<lwe::Integer> result(num_rows, 0);
  for (int j = 0; j < vec.size(); ++j) {
    int i = 0;
    for (int block_idx = 0; block_idx < num_blocks; ++block_idx) {
      BlockType block = matrix[j][block_idx];
      const DBIntegerType* block_as_values =
          reinterpret_cast<const DBIntegerType*>(&block);
      // NOTE: This reinterpret cast is very unstable. It only really words for u8. For u16, it will not work, but we know that our database values are at most 8-bits, so we can use this for now.
      for (int block_pos = 0; block_pos < num_values_per_block && i < num_rows;
           ++block_pos, ++i) {
            result[i] +=
              (static_cast<lwe::Integer>(block_as_values[block_pos])) * vec[j];
      }
    }
  }
  return result;
}


// Only instantiate the 8-bit and 16-bit versions, which are the choices of
// LWE plaintext integer types we support.
HWY_EXPORT_T(InnerProductHwy8, InnerProductHwy<uint8_t>);
HWY_EXPORT_T(InnerProductHwy16, InnerProductHwy<uint16_t>);

HWY_EXPORT_T(InnerProductHwy8_64, InnerProductHwy64<uint8_t>);

template <>
absl::StatusOr<std::vector<lwe::Integer>> InnerProduct<uint8_t>(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec, bool secondary) {
  return HWY_DYNAMIC_DISPATCH_T(InnerProductHwy8)(matrix, vec, secondary);
}

template <>
absl::StatusOr<std::vector<lwe::Integer>> InnerProduct<uint16_t>(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer> vec, bool secondary) {
  return HWY_DYNAMIC_DISPATCH_T(InnerProductHwy16)(matrix, vec, secondary);
}

template <>
absl::StatusOr<std::vector<lwe::Integer64>> InnerProduct64<uint8_t>(
    absl::Span<const BlockVector> matrix, absl::Span<const lwe::Integer64> vec, bool secondary) {
  return HWY_DYNAMIC_DISPATCH_T(InnerProductHwy8_64)(matrix, vec, secondary);
}

}  // namespace hintless_pir::hintless_simplepir::internal
#endif  // HWY_ONCE || HWY_IDE
