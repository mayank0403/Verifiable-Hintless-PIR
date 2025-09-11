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

#include "hintless_simplepir/database_hwy.h"

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/memory/memory.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "hintless_simplepir/inner_product_hwy.h"
#include "hintless_simplepir/parameters.h"
#include "hintless_simplepir/utils.h"
#include "lwe/types.h"
#include "shell_encryption/status_macros.h"

namespace hintless_pir {
namespace hintless_simplepir {
namespace {

static inline Database::RawMatrix CreateZeroRawMatrix(size_t num_rows,
                                                      size_t num_cols) {
  //size_t num_values_per_block =
  //    sizeof(internal::BlockType) / sizeof(lwe::PlainInteger);
  size_t num_values_per_block =
      sizeof(internal::BlockType) / sizeof(lwe::DBInteger);
  size_t num_blocks_per_col = DivAndRoundUp(num_rows, num_values_per_block);
  Database::RawMatrix matrix(num_cols);
  for (int i = 0; i < num_cols; ++i) {
    matrix[i].resize(num_blocks_per_col, 0);
  }
  return matrix;
}

static inline Database::RawMatrix CreateRandomRawMatrix(size_t num_rows,
                                                        size_t num_cols,
                                                        size_t plain_bits) {
  //size_t num_values_per_block =
  //    sizeof(internal::BlockType) / sizeof(lwe::PlainInteger);
  size_t num_values_per_block =
      sizeof(internal::BlockType) / sizeof(lwe::DBInteger);
  size_t num_blocks_per_col = DivAndRoundUp(num_rows, num_values_per_block);
  lwe::Integer mask = (lwe::Integer{1} << plain_bits) - 1;
  Database::RawMatrix matrix(num_cols);
  for (int i = 0; i < num_cols; ++i) {
    matrix[i].resize(num_blocks_per_col, 0);
    for (int j = 0; j < num_blocks_per_col; ++j) {
      //for (int k = 0, b = 0; k < num_values_per_block;
      //     ++k, b += 8 * sizeof(lwe::PlainInteger)) {
      for (int k = 0, b = 0; k < num_values_per_block;
           ++k, b += 8 * sizeof(lwe::DBInteger)) {
        lwe::Integer r = std::rand();
        matrix[i][j] |= static_cast<internal::BlockType>(r & mask) << b;
      }
    }
  }
  return matrix;
}

// This function will be quite slow for large matrices because it requires a transpose which is makes many random-ish accesses.
static inline std::pair<Database::RawMatrix, Database::RawMatrix> CreateRandomRawMatrixWithTransposeFake(size_t num_rows,
                                                        size_t num_cols,
                                                        size_t plain_bits) {
  //size_t num_values_per_block =
  //    sizeof(internal::BlockType) / sizeof(lwe::PlainInteger);
  size_t num_values_per_block =
      sizeof(internal::BlockType) / sizeof(lwe::DBInteger);
  size_t num_blocks_per_col = DivAndRoundUp(num_rows, num_values_per_block);
  size_t num_blocks_per_row = DivAndRoundUp(num_cols, num_values_per_block);
  lwe::Integer mask = (lwe::Integer{1} << plain_bits) - 1;
  Database::RawMatrix matrix(num_cols);
  Database::RawMatrix matrix_transpose(num_rows);

  // Column form
  //Database::LweMatrix uncompr_matrix(num_cols);
  //Database::LweMatrix uncompr_matrix_transpose(num_rows);

  // Initialize matrix
  for (int i = 0; i < num_cols; ++i) {
    matrix[i].resize(num_blocks_per_col, 0);
  }

  // Initialize matrix_transpose
  for (int i = 0; i < num_rows; i++) {
    matrix_transpose[i].resize(num_blocks_per_row, 0);
  }

  // Tracking current location in uncompressed matrix
  size_t cur_row_uncmpr = 0;
  size_t cur_col_uncmpr = 0;

  for (int i = 0; i < num_cols; ++i) {
    cur_row_uncmpr = 0;
    cur_col_uncmpr = i;
    for (int j = 0; j < num_blocks_per_col; ++j) {
      //for (int k = 0, b = 0; k < num_values_per_block;
      //     ++k, b += 8 * sizeof(lwe::PlainInteger)) {
      for (int k = 0, b = 0; k < num_values_per_block;
           ++k, b += 8 * sizeof(lwe::DBInteger)) {
        lwe::Integer r = std::rand();
        // Write to matrix
        matrix[i][j] |= static_cast<internal::BlockType>(r & mask) << b;
        cur_row_uncmpr += 1;
        if (cur_row_uncmpr >= num_rows) {
          break;
        }
      }
    }
  }

    // Handles matrix_transpose only
  for (int i = 0; i < num_rows; ++i) {
    cur_row_uncmpr = i;
    cur_col_uncmpr = 0;
    for (int j = 0; j < num_blocks_per_row; ++j) {
      //for (int k = 0, b = 0; k < num_values_per_block;
      //     ++k, b += 8 * sizeof(lwe::PlainInteger)) {
      for (int k = 0, b = 0; k < num_values_per_block;
           ++k, b += 8 * sizeof(lwe::DBInteger)) {
        lwe::Integer r = std::rand();
        // Write to matrix
        matrix_transpose[i][j] |= static_cast<internal::BlockType>(r & mask) << b;
        cur_col_uncmpr += 1;
        if (cur_col_uncmpr >= num_cols) {
          break;
        }
      }
    }
  }

  return make_pair(matrix, matrix_transpose);
}

// This function will be quite slow for large matrices because it requires a transpose which is makes many random-ish accesses.
static inline std::pair<Database::RawMatrix, Database::RawMatrix> CreateRandomRawMatrixWithTranspose(size_t num_rows,
                                                        size_t num_cols,
                                                        size_t plain_bits) {
  //size_t num_values_per_block =
  //    sizeof(internal::BlockType) / sizeof(lwe::PlainInteger);
  size_t num_values_per_block =
      sizeof(internal::BlockType) / sizeof(lwe::DBInteger);
  size_t num_blocks_per_col = DivAndRoundUp(num_rows, num_values_per_block);
  size_t num_blocks_per_row = DivAndRoundUp(num_cols, num_values_per_block);
  lwe::Integer mask = (lwe::Integer{1} << plain_bits) - 1;
  Database::RawMatrix matrix(num_cols);
  Database::RawMatrix matrix_transpose(num_rows);

  // Column form
  //Database::LweMatrix uncompr_matrix(num_cols);
  //Database::LweMatrix uncompr_matrix_transpose(num_rows);

  // Initialize matrix
  for (int i = 0; i < num_cols; ++i) {
    matrix[i].resize(num_blocks_per_col, 0);
  }

  // Initialize matrix_transpose
  for (int i = 0; i < num_rows; i++) {
    matrix_transpose[i].resize(num_blocks_per_row, 0);
  }

  // Tracking current location in uncompressed matrix
  size_t cur_row_uncmpr = 0;
  size_t cur_col_uncmpr = 0;

  for (int i = 0; i < num_cols; ++i) {
    cur_row_uncmpr = 0;
    cur_col_uncmpr = i;
    size_t transpose_block_id = cur_col_uncmpr / num_values_per_block;
    size_t transpose_block_offset = cur_col_uncmpr % num_values_per_block;
    for (int j = 0; j < num_blocks_per_col; ++j) {
      //for (int k = 0, b = 0; k < num_values_per_block;
      //     ++k, b += 8 * sizeof(lwe::PlainInteger)) {
      for (int k = 0, b = 0; k < num_values_per_block;
           ++k, b += 8 * sizeof(lwe::DBInteger)) {
        lwe::Integer r = std::rand();
        // Write to matrix
        matrix[i][j] |= static_cast<internal::BlockType>(r & mask) << b;
        // Write to matrix_transpose
        matrix_transpose[cur_row_uncmpr][transpose_block_id] |= static_cast<internal::BlockType>(r & mask) << (8 * sizeof(lwe::DBInteger) * transpose_block_offset);
        cur_row_uncmpr += 1;
        if (cur_row_uncmpr >= num_rows) {
          break;
        }
      }
    }
  }
  return make_pair(matrix, matrix_transpose);
}

// This function is quite fast. It just put static values in the matrix. At location (i, j) it puts i + 3*j.
static inline std::pair<Database::RawMatrix, Database::RawMatrix> CreateFastRawMatrixWithTranspose(size_t num_rows,
                                                        size_t num_cols,
                                                        size_t plain_bits) {
  //size_t num_values_per_block =
  //    sizeof(internal::BlockType) / sizeof(lwe::PlainInteger);
  size_t num_values_per_block =
      sizeof(internal::BlockType) / sizeof(lwe::DBInteger);
  size_t num_blocks_per_col = DivAndRoundUp(num_rows, num_values_per_block);
  size_t num_blocks_per_row = DivAndRoundUp(num_cols, num_values_per_block);
  lwe::Integer mask = (lwe::Integer{1} << plain_bits) - 1;
  Database::RawMatrix matrix(num_cols);
  Database::RawMatrix matrix_transpose(num_rows);

  // Column form
  //Database::LweMatrix uncompr_matrix(num_cols);
  //Database::LweMatrix uncompr_matrix_transpose(num_rows);

  // Initialize matrix
  for (int i = 0; i < num_cols; ++i) {
    matrix[i].resize(num_blocks_per_col, 0);
  }

  // Initialize matrix_transpose
  for (int i = 0; i < num_rows; i++) {
    matrix_transpose[i].resize(num_blocks_per_row, 0);
  }

  // Tracking current location in uncompressed matrix
  size_t cur_row_uncmpr = 0;
  size_t cur_col_uncmpr = 0;

  // Handles matrix only
  for (int i = 0; i < num_cols; ++i) {
    cur_row_uncmpr = 0;
    cur_col_uncmpr = i;
    for (int j = 0; j < num_blocks_per_col; ++j) {
      //for (int k = 0, b = 0; k < num_values_per_block;
      //     ++k, b += 8 * sizeof(lwe::PlainInteger)) {
      for (int k = 0, b = 0; k < num_values_per_block;
           ++k, b += 8 * sizeof(lwe::DBInteger)) {
        lwe::Integer r = cur_row_uncmpr + 3 * cur_col_uncmpr;
        // Write to matrix
        matrix[i][j] |= static_cast<internal::BlockType>(r & mask) << b;
        cur_row_uncmpr += 1;
        if (cur_row_uncmpr >= num_rows) {
          break;
        }
      }
    }
  }

    // Handles matrix_transpose only
  for (int i = 0; i < num_rows; ++i) {
    cur_row_uncmpr = i;
    cur_col_uncmpr = 0;
    for (int j = 0; j < num_blocks_per_row; ++j) {
      //for (int k = 0, b = 0; k < num_values_per_block;
      //     ++k, b += 8 * sizeof(lwe::PlainInteger)) {
      for (int k = 0, b = 0; k < num_values_per_block;
           ++k, b += 8 * sizeof(lwe::DBInteger)) {
        lwe::Integer r = cur_row_uncmpr + 3 * cur_col_uncmpr;
        // Write to matrix
        matrix_transpose[i][j] |= static_cast<internal::BlockType>(r & mask) << b;
        cur_col_uncmpr += 1;
        if (cur_col_uncmpr >= num_cols) {
          break;
        }
      }
    }
  }

  return make_pair(matrix, matrix_transpose);
}


static inline Database::LweMatrix CreateZeroMatrix(size_t num_rows,
                                                   size_t num_cols) {
  Database::LweMatrix matrix(num_rows);
  for (int i = 0; i < num_rows; ++i) {
    matrix[i].resize(num_cols, 0);
  }
  return matrix;
}

static inline Database::LweMatrix64 CreateZeroMatrix64(size_t num_rows,
                                                   size_t num_cols) {
  Database::LweMatrix64 matrix(num_rows);
  for (int i = 0; i < num_rows; ++i) {
    matrix[i].resize(num_cols, 0);
  }
  return matrix;
}

// Assume both `plain_matrix` and `lwe_matrix` are stored by columns.
static inline absl::StatusOr<Database::LweMatrix> MatrixProduct(
    const Database::RawMatrix& plain_matrix,
    const Database::LweMatrix& lwe_matrix, size_t num_rows, bool secondary = false) {
  Database::LweMatrix cols;
  cols.reserve(lwe_matrix.size());
  for (int i = 0; i < lwe_matrix.size(); ++i) {
    /*RLWE_ASSIGN_OR_RETURN(
        Database::LweVector col,
        internal::InnerProduct<lwe::PlainInteger>(plain_matrix, lwe_matrix[i], solinas, (i<=0)));*/
    RLWE_ASSIGN_OR_RETURN(
        Database::LweVector col,
        internal::InnerProduct<lwe::DBInteger>(plain_matrix, lwe_matrix[i], secondary));
    cols.push_back(col);
  }
  // return `matrix` organized by rows.
  Database::LweMatrix matrix(num_rows);
  for (int i = 0; i < num_rows; ++i) {
    matrix[i].resize(lwe_matrix.size());
    for (int j = 0; j < lwe_matrix.size(); ++j) {
      matrix[i][j] = cols[j][i];
    }
  }
  return matrix;
}

// Assume both `plain_matrix` and `lwe_matrix` are stored by columns.
static inline absl::StatusOr<Database::LweMatrix64> MatrixProduct64(
    const Database::RawMatrix& plain_matrix,
    const Database::LweMatrix64& lwe_matrix, size_t num_rows, bool secondary = false) {
  Database::LweMatrix64 cols;
  cols.reserve(lwe_matrix.size());
  for (int i = 0; i < lwe_matrix.size(); ++i) {
    /*RLWE_ASSIGN_OR_RETURN(
        Database::LweVector64 col,
        internal::InnerProduct64<lwe::PlainInteger>(plain_matrix, lwe_matrix[i], solinas, (i<=0)));*/
    RLWE_ASSIGN_OR_RETURN(
        Database::LweVector64 col,
        internal::InnerProduct64<lwe::DBInteger>(plain_matrix, lwe_matrix[i], secondary));
    cols.push_back(col);
  }
  // return `matrix` organized by rows.
  Database::LweMatrix64 matrix(num_rows);
  for (int i = 0; i < num_rows; ++i) {
    matrix[i].resize(lwe_matrix.size());
    for (int j = 0; j < lwe_matrix.size(); ++j) {
      matrix[i][j] = cols[j][i];
    }
  }
  return matrix;
}

}  // namespace

absl::StatusOr<std::unique_ptr<Database>> Database::Create(
    const Parameters& parameters) {
  // Initialize the data and the hint matrices for all shards.
  int num_shards = DivAndRoundUp(parameters.db_record_bit_size,
                                 parameters.lwe_plaintext_bit_size);
  std::vector<RawMatrix> data_matrices(num_shards);
  std::vector<RawMatrix> data_matrices_transpose(num_shards);
  std::vector<LweMatrix> hint_matrices(num_shards);
  for (int i = 0; i < num_shards; ++i) {
    data_matrices[i] =
        CreateZeroRawMatrix(parameters.db_rows, parameters.db_cols);
    data_matrices_transpose[i] =
        CreateZeroRawMatrix(parameters.db_cols, parameters.db_rows);
    hint_matrices[i] =
        CreateZeroMatrix(parameters.db_rows, parameters.lwe_secret_dim);
  }
  return absl::WrapUnique(
      new Database(parameters, /*lwe_query_pad=*/nullptr, /*num_records=*/0,
                   std::move(data_matrices), std::move(data_matrices_transpose), std::move(hint_matrices)));
}

absl::StatusOr<std::unique_ptr<Database>> Database::CreateRandomFake(
    const Parameters& parameters) {
  // Initialize the data and the hint matrices for all shards.
  int num_shards = DivAndRoundUp(parameters.db_record_bit_size,
                                 parameters.lwe_plaintext_bit_size);
  int db_record_bits = parameters.lwe_plaintext_bit_size;
  if (db_record_bits > parameters.db_record_bit_size) {
    db_record_bits = parameters.db_record_bit_size;
  }
  //std::cout << "db_record_bits: " << db_record_bits << std::endl;
  //std::cout << "num_shards: " << num_shards << std::endl;
  std::vector<RawMatrix> data_matrices(num_shards);
  std::vector<RawMatrix> data_matrices_transpose(num_shards);
  std::vector<LweMatrix> hint_matrices(num_shards);
  for (int i = 0; i < num_shards; ++i) {
    auto mats =
        CreateRandomRawMatrixWithTransposeFake(parameters.db_rows, parameters.db_cols,
                              db_record_bits);
    data_matrices[i] = mats.first;
    data_matrices_transpose[i] = mats.second;
    hint_matrices[i] =
        CreateZeroMatrix(parameters.db_rows, parameters.lwe_secret_dim);
  }
  int64_t num_records = parameters.db_rows * parameters.db_cols;
  return absl::WrapUnique(new Database(parameters, /*lwe_query_pad=*/nullptr,
                                       num_records, std::move(data_matrices), std::move(data_matrices_transpose),
                                       std::move(hint_matrices)));
}

absl::StatusOr<std::unique_ptr<Database>> Database::CreateRandom(
    const Parameters& parameters) {
  // Initialize the data and the hint matrices for all shards.
  int num_shards = DivAndRoundUp(parameters.db_record_bit_size,
                                 parameters.lwe_plaintext_bit_size);
  int db_record_bits = parameters.lwe_plaintext_bit_size;
  if (db_record_bits > parameters.db_record_bit_size) {
    db_record_bits = parameters.db_record_bit_size;
  }
  //std::cout << "db_record_bits: " << db_record_bits << std::endl;
  //std::cout << "num_shards: " << num_shards << std::endl;
  std::vector<RawMatrix> data_matrices(num_shards);
  std::vector<RawMatrix> data_matrices_transpose(num_shards);
  std::vector<LweMatrix> hint_matrices(num_shards);
  for (int i = 0; i < num_shards; ++i) {
    auto mats =
        CreateRandomRawMatrixWithTranspose(parameters.db_rows, parameters.db_cols,
                              db_record_bits);
    data_matrices[i] = mats.first;
    data_matrices_transpose[i] = mats.second;
    hint_matrices[i] =
        CreateZeroMatrix(parameters.db_rows, parameters.lwe_secret_dim);
  }
  int64_t num_records = parameters.db_rows * parameters.db_cols;
  return absl::WrapUnique(new Database(parameters, /*lwe_query_pad=*/nullptr,
                                       num_records, std::move(data_matrices), std::move(data_matrices_transpose),
                                       std::move(hint_matrices)));
}

absl::StatusOr<std::unique_ptr<Database>> Database::CreateFast(
    const Parameters& parameters) {
  // Initialize the data and the hint matrices for all shards.
  int num_shards = DivAndRoundUp(parameters.db_record_bit_size,
                                 parameters.lwe_plaintext_bit_size);
  int db_record_bits = parameters.lwe_plaintext_bit_size;
  if (db_record_bits > parameters.db_record_bit_size) {
    db_record_bits = parameters.db_record_bit_size;
  }
  //std::cout << "db_record_bits: " << db_record_bits << std::endl;
  //std::cout << "num_shards: " << num_shards << std::endl;
  std::vector<RawMatrix> data_matrices(num_shards);
  std::vector<RawMatrix> data_matrices_transpose(num_shards);
  std::vector<LweMatrix> hint_matrices(num_shards);
  for (int i = 0; i < num_shards; ++i) {
    auto mats =
        CreateFastRawMatrixWithTranspose(parameters.db_rows, parameters.db_cols,
                              db_record_bits);
    data_matrices[i] = mats.first;
    data_matrices_transpose[i] = mats.second;
    hint_matrices[i] =
        CreateZeroMatrix(parameters.db_rows, parameters.lwe_secret_dim);
  }
  int64_t num_records = parameters.db_rows * parameters.db_cols;
  return absl::WrapUnique(new Database(parameters, /*lwe_query_pad=*/nullptr,
                                       num_records, std::move(data_matrices), std::move(data_matrices_transpose),
                                       std::move(hint_matrices)));
}

absl::Status Database::UpdateLweQueryPad(const lwe::Matrix* lwe_query_pad) {
  if (lwe_query_pad == nullptr) {
    return absl::InvalidArgumentError("`lwe_query_pad` must not be null.");
  }
  lwe_query_pad_ = lwe_query_pad;
  return absl::OkStatus();
}

absl::Status Database::Append(absl::string_view record) {
  if (record.size() * 8 >= params_.db_record_bit_size + 8 ||
      record.size() * 8 < params_.db_record_bit_size) {
    return absl::InvalidArgumentError("`record` has incorrect size.");
  }
  if (num_records_ >= params_.db_rows * params_.db_cols) {
    return absl::InvalidArgumentError("Database is full.");
  }
  int64_t row_idx, col_idx;
  std::tie(row_idx, col_idx) = MatrixCoordinate(num_records_);
  //int64_t num_values_per_block = sizeof(BlockType) / sizeof(lwe::PlainInteger);
  int64_t num_values_per_block = sizeof(BlockType) / sizeof(lwe::DBInteger);
  int64_t block_idx = row_idx / num_values_per_block;
  int64_t block_pos = row_idx % num_values_per_block;
  //int64_t base_bits = block_pos * 8 * sizeof(lwe::PlainInteger);
  int64_t base_bits = block_pos * 8 * sizeof(lwe::DBInteger);

  num_records_++;
  std::vector<lwe::Integer> values = SplitRecord(record, params_);
  for (int i = 0; i < values.size(); ++i) {
    BlockType block = static_cast<BlockType>(values[i]) << base_bits;
    data_matrices_[i][col_idx][block_idx] |= block;
  }
  return absl::OkStatus();
}

absl::Status Database::UpdateHints() {
  if (lwe_query_pad_ == nullptr) {
    return absl::FailedPreconditionError("LWE query pad not set.");
  }
  for (int i = 0; i < data_matrices_.size(); ++i) {
    LweMatrix lwe_matrix = ImportLweMatrix(*lwe_query_pad_);
    /*
    std::cout << "Cols" << lwe_matrix.size() << " Rows" << lwe_matrix[0].size() << std::endl;
    std::cout << "Hintless A 0th column ";
    for (int j = 0; j < lwe_matrix[0].size(); ++j) {
      std::cout << " " << lwe_matrix[0][j];
    }
    std::cout << std::endl;
    */
    
    RLWE_ASSIGN_OR_RETURN(
        hint_matrices_[i],
        MatrixProduct(data_matrices_[i], lwe_matrix, params_.db_rows));
  }
  return absl::OkStatus();
}

absl::Status Database::UpdateHintsFake() {
  std::cout << "Generating FAKE hint" << std::endl;
  if (lwe_query_pad_ == nullptr) {
    return absl::FailedPreconditionError("LWE query pad not set.");
  }
  int rows = params_.db_rows;
  int cols = params_.lwe_secret_dim;
  for (int i = 0; i < data_matrices_.size(); ++i) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        hint_matrices_[i][r][c] = static_cast<lwe::Integer>(std::rand()) * static_cast<lwe::Integer>(std::rand()); // To ensure result is at 64-bits on all machines
      }
    }
  }
  return absl::OkStatus();
}

// Assume both `plain_matrix` and `lwe_matrix` are stored by columns.
absl::StatusOr<std::vector<Database::LweMatrix>> Database::TransposeMatrixProduct(
    const Database::LweMatrix& lwe_matrix, bool secondary) const {
  size_t shards = data_matrices_transpose_.size();
  std::vector<Database::LweMatrix> res(shards);
  for (int i = 0; i < shards; ++i) {
    RLWE_ASSIGN_OR_RETURN(
        res[i],
        MatrixProduct(data_matrices_transpose_[i], lwe_matrix, params_.db_cols, secondary));
  }
  return res;
}

// Assume both `plain_matrix` and `lwe_matrix` are stored by columns.
absl::StatusOr<std::vector<Database::LweMatrix64>> Database::TransposeMatrixProduct64(
    const Database::LweMatrix64& lwe_matrix, bool secondary) const {
  size_t shards = data_matrices_transpose_.size();
  std::vector<Database::LweMatrix64> res(shards);
  for (int i = 0; i < shards; ++i) {
    RLWE_ASSIGN_OR_RETURN(
        res[i],
        MatrixProduct64(data_matrices_transpose_[i], lwe_matrix, params_.db_cols, secondary));
  }
  return res;
}

// Assume both `plain_matrix` and `lwe_matrix` are stored by columns.
absl::StatusOr<std::vector<Database::LweMatrix>> Database::DirectMatrixProduct(
    const Database::LweMatrix& lwe_matrix, bool secondary) const {
  size_t shards = data_matrices_.size();
  std::vector<Database::LweMatrix> res(shards);
  for (int i = 0; i < shards; ++i) {
    RLWE_ASSIGN_OR_RETURN(
        res[i],
        MatrixProduct(data_matrices_[i], lwe_matrix, params_.db_rows, secondary));
  }
  return res;
}

std::vector<Database::LweMatrix> Database::TransposeMatrixProductFake() const {
  std::cout << "Fake Hint generation..." << std::endl;
  size_t shards = data_matrices_transpose_.size();
  std::vector<Database::LweMatrix> res(shards);
  int rows = params_.db_cols;
  int cols = params_.lwe_secret_dim;
  for (int i = 0; i < shards; ++i) {
    res[i] = CreateZeroMatrix(rows, cols);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        res[i][r][c] = static_cast<lwe::Integer>(std::rand()) * static_cast<lwe::Integer>(std::rand()); // To ensure result is at 64-bits on all machines
      }
    }
  }
  return res;
}

std::vector<Database::LweMatrix64> Database::TransposeMatrixProductFake64() const {
  std::cout << "Fake (64-bit size) Hint generation..." << std::endl;
  size_t shards = data_matrices_transpose_.size();
  std::vector<Database::LweMatrix64> res(shards);
  int rows = params_.db_cols;
  int cols = params_.offline_lwe_secret_dim;
  for (int i = 0; i < shards; ++i) {
    res[i] = CreateZeroMatrix64(rows, cols);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        res[i][r][c] = static_cast<lwe::Integer64>(std::rand()) * static_cast<lwe::Integer64>(std::rand()); // To ensure result is at 64-bits on all machines
      }
    }
  }
  return res;
}

std::vector<Database::LweMatrix> Database::DirectMatrixProductFake() const {
  std::cout << "Fake Hint generation..." << std::endl;
  size_t shards = data_matrices_.size();
  std::vector<Database::LweMatrix> res(shards);
  int rows = params_.db_rows;
  int cols = params_.lwe_secret_dim;
  for (int i = 0; i < shards; ++i) {
    res[i] = CreateZeroMatrix(rows, cols);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        res[i][r][c] = static_cast<lwe::Integer>(std::rand()) * static_cast<lwe::Integer>(std::rand()); // To ensure result is at 64-bits on all machines
      }
    }
  }
  return res;
}

// Assume both `plain_matrix` is stored by columns.
absl::StatusOr<std::vector<std::vector<Database::LweVector>>> Database::TransposeMatrixCiphertextProducts(
    const std::vector<Database::LweVector>& cts) const {
  size_t shards = data_matrices_transpose_.size();
  std::vector<std::vector<Database::LweVector>> res(shards);
  for (int i = 0; i < shards; ++i) {
    res[i].reserve(cts.size());
    for (int j = 0; j < cts.size(); j++){
      /*RLWE_ASSIGN_OR_RETURN(
          Database::LweVector col,
          internal::InnerProduct<lwe::PlainInteger>(plain_matrix, lwe_matrix[i], solinas, (i<=0)));*/
      RLWE_ASSIGN_OR_RETURN(
          Database::LweVector col,
          internal::InnerProduct<lwe::DBInteger>(data_matrices_transpose_[i], cts[j]));
      res[i].push_back(col);
    }
  }
  return res;
}

// Assume both `plain_matrix` is stored by columns.
absl::StatusOr<std::vector<std::vector<Database::LweVector64>>> Database::TransposeMatrixCiphertextProducts64(
    const std::vector<Database::LweVector64>& cts) const {
  size_t shards = data_matrices_transpose_.size();
  std::vector<std::vector<Database::LweVector64>> res(shards);
  for (int i = 0; i < shards; ++i) {
    res[i].reserve(cts.size());
    for (int j = 0; j < cts.size(); j++){
      /*RLWE_ASSIGN_OR_RETURN(
          Database::LweVector64 col,
          internal::InnerProduct64<lwe::PlainInteger>(plain_matrix, lwe_matrix[i], solinas, (i<=0)));*/
      RLWE_ASSIGN_OR_RETURN(
          Database::LweVector64 col,
          internal::InnerProduct64<lwe::DBInteger>(data_matrices_transpose_[i], cts[j]));
      res[i].push_back(col);
    }
  }
  return res;
}

absl::StatusOr<std::vector<Database::LweVector>> Database::InnerProductWith(
    const LweVector& query) const {
  std::vector<LweVector> results;
  results.reserve(data_matrices_.size());
  for (auto const& matrix : data_matrices_) {
    /*RLWE_ASSIGN_OR_RETURN(
        LweVector result,
        internal::InnerProduct<lwe::PlainInteger>(matrix, query, solinas));*/
    RLWE_ASSIGN_OR_RETURN(
        LweVector result,
        internal::InnerProduct<lwe::DBInteger>(matrix, query));
    result.resize(params_.db_rows);
    results.push_back(result);
  }
  return results;
}

/*
absl::StatusOr<std::vector<Database::LweVector>> Database::InnerProductWithFaster32by8(
    const LweVector& query) const {
  std::cout << "[~~ OPTIMIZATION ~~]: Faster 32+16 bit inner product in online phase" << std::endl;
  std::vector<LweVector> results;
  results.reserve(data_matrices_.size());
  std::vector<lwe::Integer> query_upper; // upper 24 bits
  std::vector<uint16_t> query_lower; // lower 8 bits
  lwe::Integer mask = (lwe::Integer{1} << 8) - 1;
  for (int i = 0; i < query.size(); i++) {
    query_upper.push_back(query[i] >> 8);
    query_lower.push_back(static_cast<uint16_t>(query[i] & mask));
  }
  for (auto const& matrix : data_matrices_) {
    RLWE_ASSIGN_OR_RETURN(
        LweVector result,
        internal::InnerProductFast32by8<lwe::DBInteger>(matrix, query_upper, query_lower));
    result.resize(params_.db_rows);
    results.push_back(result);
  }
  return results;
}
*/

absl::StatusOr<std::string> Database::Record(int64_t index) const {
  if (index < 0 || index >= num_records_) {
    return absl::InvalidArgumentError("`index` is out of range.");
  }
  int64_t row_idx, col_idx;
  std::tie(row_idx, col_idx) = MatrixCoordinate(index);
  //int64_t num_values_per_block = sizeof(BlockType) / sizeof(lwe::PlainInteger);
  int64_t num_values_per_block = sizeof(BlockType) / sizeof(lwe::DBInteger);
  int64_t block_idx = row_idx / num_values_per_block;
  int64_t block_pos = row_idx % num_values_per_block;
  //int64_t base_bits = block_pos * 8 * sizeof(lwe::PlainInteger);
  int64_t base_bits = block_pos * 8 * sizeof(lwe::DBInteger);

  BlockType mask = (BlockType{1} << params_.lwe_plaintext_bit_size) - 1;
  std::vector<lwe::Integer> values;
  values.reserve(data_matrices_.size());
  for (auto const& data_matrix : data_matrices_) {
    BlockType block = data_matrix[col_idx][block_idx] >> base_bits;
    values.push_back(static_cast<lwe::Integer>(block & mask));
  }
  return ReconstructRecord(values, params_);
}

absl::StatusOr<std::string> Database::RecordWithStacking(int64_t index) const {
  if (index < 0 || index >= num_records_) {
    return absl::InvalidArgumentError("`index` is out of range.");
  }
  int64_t stack_size = params_.db_stack_cells;
  int64_t row_idx, col_idx, row_idx_end;
  std::tie(row_idx, col_idx) = MatrixCoordinate(index);
  std::cout << "Row id " << row_idx << " Col id " << col_idx << " stack size " << stack_size << std::endl;
  row_idx_end = row_idx + stack_size;
  //int64_t num_values_per_block = sizeof(BlockType) / sizeof(lwe::PlainInteger);
  int64_t num_values_per_block = sizeof(BlockType) / sizeof(lwe::DBInteger);

  std::vector<std::tuple<int64_t, int64_t, int64_t>> block_locs;
  for (int64_t i = 0; i < stack_size; i++) {
    int64_t block_idx = (row_idx + i) / num_values_per_block;
    int64_t block_pos = (row_idx + i) % num_values_per_block;
    int64_t base_bits = block_pos * 8 * sizeof(lwe::DBInteger);
    block_locs.push_back(std::make_tuple(block_idx, block_pos, base_bits));
  }
  
  //int64_t block_idx_end = row_idx_end / num_values_per_block;
  //int64_t block_pos_end = row_idx_end % num_values_per_block;
  //int64_t base_bits = block_pos * 8 * sizeof(lwe::PlainInteger);
  //int64_t base_bits = block_pos * 8 * sizeof(lwe::DBInteger);

  BlockType mask = (BlockType{1} << params_.lwe_plaintext_bit_size) - 1;
  std::vector<lwe::Integer> values;
  values.reserve(data_matrices_.size() * stack_size);
  for (auto const& data_matrix : data_matrices_) {
    for (int64_t i = 0; i < stack_size; i++) {
      //BlockType block = data_matrix[col_idx][block_idx] >> base_bits;
      BlockType block = data_matrix[col_idx][std::get<0>(block_locs[i])] >> std::get<2>(block_locs[i]);
      values.push_back(static_cast<lwe::Integer>(block & mask));
      //std::cout << i << " " << static_cast<lwe::Integer>(block & mask) << std::endl;
    }
  }
  return ReconstructRecordWithStacking(values, params_);
}

Database::LweMatrix ImportLweMatrix(const lwe::Matrix& matrix) {
  // `results` organized by columns.
  Database::LweMatrix results(matrix.cols());
  for (int64_t j = 0; j < matrix.cols(); ++j) {
    results[j].resize(matrix.rows(), 0);
    for (int64_t i = 0; i < matrix.rows(); ++i) {
      results[j][i] = matrix(i, j);
    }
  }
  return results;
}

lwe::Matrix ExportLweMatrix(const Database::LweMatrix& matrix) {
  // Assume `matrix` organized by columns.
  int64_t num_cols = matrix.size();
  int64_t num_rows = matrix[0].size();
  lwe::Matrix results = lwe::Matrix::Zero(num_rows, num_cols);
  for (int64_t j = 0; j < num_cols; ++j) {
    for (int64_t i = 0; i < num_rows; ++i) {
      results(i, j) = matrix[j][i];
    }
  }
  return results;
}

lwe::Matrix ExportRawMatrix(const Database::RawMatrix& matrix, size_t num_rows,
                            size_t num_bits_per_value) {
  // Assume `matrix` organized by columns.
  int64_t num_cols = matrix.size();
  //int64_t num_values_per_block =
  //    Database::kBlockBits / sizeof(lwe::PlainInteger);
  int64_t num_values_per_block =
      Database::kBlockBits / sizeof(lwe::DBInteger);
  lwe::Integer mask = (lwe::Integer{1} << num_bits_per_value) - 1;
  lwe::Matrix results = lwe::Matrix::Zero(num_rows, num_cols);
  for (int64_t col_idx = 0; col_idx < num_cols; ++col_idx) {
    for (int64_t row_idx = 0; row_idx < num_rows; ++row_idx) {
      int64_t block_idx = row_idx / num_values_per_block;
      int64_t block_pos = row_idx % num_values_per_block;
      //int64_t base_bits = block_pos * 8 * sizeof(lwe::PlainInteger);
      int64_t base_bits = block_pos * 8 * sizeof(lwe::DBInteger);
      auto raw =
          static_cast<lwe::Integer>(matrix[col_idx][block_idx] >> base_bits);
      results(row_idx, col_idx) = raw & mask;
    }
  }
  return results;
}

}  // namespace hintless_simplepir
}  // namespace hintless_pir
