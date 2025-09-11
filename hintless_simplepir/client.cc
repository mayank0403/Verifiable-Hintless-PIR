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

#include "hintless_simplepir/client.h"

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Eigen/Core"
#include "absl/memory/memory.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "hintless_simplepir/parameters.h"
#include "hintless_simplepir/serialization.pb.h"
#include "hintless_simplepir/utils.h"
#include "lwe/encode.h"
#include "lwe/lwe_symmetric_encryption.h"
#include "lwe/types.h"
#include "shell_encryption/int256.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/crt_interpolation.h"
#include "shell_encryption/status_macros.h"

namespace hintless_pir {
namespace hintless_simplepir {


absl::StatusOr<std::unique_ptr<Client>> Client::Create(
    const Parameters& params,
    const HintlessPirServerPublicParams& public_params) {
  if (!(params.prng_type == rlwe::PRNG_TYPE_HKDF ||
        params.prng_type == rlwe::PRNG_TYPE_CHACHA)) {
    return absl::InvalidArgumentError("Invalid PRNG type in `params`.");
  }
  // Create LinPir clients, one per plaintext modulus in `ts`.
  auto const& rlwe_params = params.linpir_params;
  int num_linpir_instances = rlwe_params.ts.size();
  if (public_params.prng_seed_linpir_ct_pads_size() != num_linpir_instances) {
    return absl::InvalidArgumentError(
        "`public_params` contains incorrect number of PRNG seeds.");
  }
  std::vector<std::unique_ptr<const RlweRnsContext>> rlwe_contexts;
  std::vector<std::unique_ptr<LinPirClient>> linpir_clients;
  rlwe_contexts.reserve(num_linpir_instances);
  linpir_clients.reserve(num_linpir_instances);
  for (int i = 0; i < num_linpir_instances; ++i) {
    RLWE_ASSIGN_OR_RETURN(
        auto rlwe_context,
        RlweRnsContext::CreateForBfvFiniteFieldEncoding(
            rlwe_params.log_n, rlwe_params.qs, /*ps=*/{}, rlwe_params.ts[i]));
    auto rlwe_context_ptr =
        std::make_unique<const RlweRnsContext>(std::move(rlwe_context));
    RLWE_ASSIGN_OR_RETURN(
        auto linpir_client,
        LinPirClient::CreateBSGS(rlwe_params, rlwe_context_ptr.get(),
                             public_params.prng_seed_linpir_ct_pads(i),
                             public_params.prng_seed_linpir_gk_pad(),
                             public_params.prng_seed_linpir_gk_pad_baby(),
                             public_params.prng_seed_linpir_gk_pad_giant()));
    rlwe_contexts.push_back(std::move(rlwe_context_ptr));
    linpir_clients.push_back(std::move(linpir_client));
  }
  std::vector<const RlwePrimeModulus*> rlwe_moduli =
      rlwe_contexts[0]->MainPrimeModuli();

  // Create another RNS context for computing CRT interpolation wrt plaintext
  // moduli. We set its plaintext modulus to 2 as a place holder.
  RLWE_ASSIGN_OR_RETURN(
      RlweRnsContext crt_context,
      RlweRnsContext::Create(rlwe_params.log_n, rlwe_params.ts, /*ps=*/{}, 2));

  return absl::WrapUnique(
      new Client(params, public_params.prng_seed_lwe_query_pad(),
                 std::move(rlwe_contexts), std::move(rlwe_moduli),
                 std::move(linpir_clients), std::move(crt_context)));
}


absl::StatusOr<HintlessPirRequest> Client::GenerateRequest(int64_t index) {
  if (index < 0 || index >= params_.db_rows * params_.db_cols) {
    return absl::InvalidArgumentError("`index` out of range.");
  }

  // Step 1. Encrypting the selection vector under LWE.
  lwe::Matrix lwe_pad;
  std::unique_ptr<rlwe::SecurePrng> lwe_enc_prng;
  std::string prng_seed_linpir_sk;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(auto pad_prng, rlwe::SingleThreadHkdfPrng::Create(
                                             prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        lwe_pad, lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim,
                                pad_prng.get()));
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadHkdfPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());

  } else {
    RLWE_ASSIGN_OR_RETURN(auto pad_prng, rlwe::SingleThreadChaChaPrng::Create(
                                             prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        lwe_pad, lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim,
                                pad_prng.get()));
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadChaChaPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
  }

  // Reduce pad modulo the prime (if prime used)
  uint64_t rows = lwe_pad.rows();
  uint64_t cols = lwe_pad.cols();
  for (uint64_t i = 0; i < rows; i++) {
    for (uint64_t j = 0; j < cols; j++) {
      lwe_pad(i, j) = fast_mod_max_64(lwe_pad(i, j));
    }
  }

  RLWE_ASSIGN_OR_RETURN(
      lwe::SymmetricLweKey lwe_secret_key,
      lwe::SymmetricLweKey::Sample(params_.lwe_secret_dim, lwe_enc_prng.get()));

  // Choosing the largest scaling factor that supports our plaintext space
  int log_scaling_factor =
      params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;

  // Plaintext is a selection vector for col_idx
  int64_t row_idx = index / params_.db_cols;
  int64_t col_idx = index % params_.db_cols;
  lwe::Vector query_vector = lwe::Vector::Zero(params_.db_cols);
  query_vector[col_idx] = 1;

  RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlace(
      query_vector, lwe_pad, log_scaling_factor, lwe_enc_prng.get()));

  // Cache the per request state.
  state_ = ClientState{.row_idx = row_idx,
                       .col_idx = col_idx,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};

  HintlessPirRequest request;
  *request.mutable_ct_query_vector() = SerializeLweCiphertext(query_vector);

  // Step 2. Encrypting the LWE secret using LinPir.
  RLWE_RETURN_IF_ERROR(
      GenerateLinPirRequestInPlace(request, lwe_secret_key.Key()));

  return request;
}

absl::StatusOr<std::pair<HintlessPirRequest, lwe::Vector>> Client::GenerateRequest(int64_t index, lwe::Vector& lwe_sk_vsp_out) {
  if (index < 0 || index >= params_.db_rows * params_.db_cols) {
    return absl::InvalidArgumentError("`index` out of range.");
  }

  // Step 1. Encrypting the selection vector under LWE.
  lwe::Matrix lwe_pad;
  std::unique_ptr<rlwe::SecurePrng> lwe_enc_prng;
  std::string prng_seed_linpir_sk;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(auto pad_prng, rlwe::SingleThreadHkdfPrng::Create(
                                             prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        lwe_pad, lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim,
                                pad_prng.get()));
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadHkdfPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());

  } else {
    RLWE_ASSIGN_OR_RETURN(auto pad_prng, rlwe::SingleThreadChaChaPrng::Create(
                                             prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        lwe_pad, lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim,
                                pad_prng.get()));
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadChaChaPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
  }

  /*
  // Reduce pad modulo the prime (if prime used)
  uint64_t rows = lwe_pad.rows();
  uint64_t cols = lwe_pad.cols();
  for (uint64_t i = 0; i < rows; i++) {
    for (uint64_t j = 0; j < cols; j++) {
      lwe_pad(i, j) = fast_mod_max_64(lwe_pad(i, j));
    }
  }
  */

  RLWE_ASSIGN_OR_RETURN(
      lwe::SymmetricLweKey lwe_secret_key,
      lwe::SymmetricLweKey::Sample(params_.lwe_secret_dim, lwe_enc_prng.get()));

  // Choosing the largest scaling factor that supports our plaintext space
  int log_scaling_factor =
      params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
  lwe::Integer scaling_factor = static_cast<lwe::Integer>(SOL_NTT_PRIME)/(1 << params_.lwe_plaintext_bit_size);
  if (scaling_factor == 0) {
    scaling_factor = (1 << log_scaling_factor);
  }

  // Plaintext is a selection vector for col_idx
  int64_t row_idx = index / params_.db_cols;
  uint64_t row_idx_end = row_idx + params_.db_stack_cells;
  int64_t col_idx = index % params_.db_cols;
  lwe::Vector query_vector = lwe::Vector::Zero(params_.db_cols);
  query_vector[col_idx] = 1;

  //RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlace(
  //    query_vector, lwe_pad, log_scaling_factor, lwe_enc_prng.get()));
  
  RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlaceExplicit(
      query_vector, lwe_pad, scaling_factor, lwe_enc_prng.get()));
  

  // Cache the per request state.
  state_ = ClientState{.row_idx = row_idx,
                       .col_idx = col_idx,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};
  state_new_ = ClientStateNew{.row_idx_begin = row_idx,
                        .row_idx_end = row_idx_end,
                       .col_idx = col_idx,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};

  HintlessPirRequest request;
  *request.mutable_ct_query_vector() = SerializeLweCiphertext(query_vector);

  double start, end;
  start = currentDateTime();
  // Step 2. Encrypting the LWE secret using LinPir.
  RLWE_RETURN_IF_ERROR(
      GenerateLinPirRequestInPlace(request, lwe_secret_key.Key()));
  end = currentDateTime();
  std::cout << "[==> TIMER  <==] Client online (LinPir) RLWE-only (rotation keys, etc.) req gen time: " << (end-start) << " ms | " << (end-start)/1000 << " sec" << std::endl;

  // Send the LWE SK to VSPIR
  for (int i = 0; i < lwe_secret_key.Key().size(); i++) {
    //std::cout << "LWE SK[" << i << "]: " << lwe_secret_key.Key()(i) << std::endl;
    lwe_sk_vsp_out(i) = lwe_secret_key.Key()(i);
    //std::cout << "LWE SK OUT[" << i << "]: " << key_out(i) << std::endl;
  }
  return std::make_pair(request, query_vector);
}

absl::StatusOr<std::pair<HintlessPirRequest, lwe::Vector>> Client::GenerateRequestGivenAs(int64_t index, lwe::Vector& As, lwe::SymmetricLweKey& s_lwe) {
  if (index < 0 || index >= params_.db_rows * params_.db_cols) {
    return absl::InvalidArgumentError("`index` out of range.");
  }

  // Step 1. Encrypting the selection vector under LWE.
  std::unique_ptr<rlwe::SecurePrng> lwe_enc_prng;
  std::string prng_seed_linpir_sk;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadHkdfPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());

  } else {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadChaChaPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
  }

  // Choosing the largest scaling factor that supports our plaintext space
  int log_scaling_factor =
      params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
  lwe::Integer scaling_factor = static_cast<lwe::Integer>(SOL_NTT_PRIME)/(1 << params_.lwe_plaintext_bit_size);
  if (scaling_factor == 0) {
    scaling_factor = (1 << log_scaling_factor);
  }

  // Plaintext is a selection vector for col_idx
  int64_t row_idx = index / params_.db_cols;
  uint64_t row_idx_end = row_idx + params_.db_stack_cells;
  int64_t col_idx = index % params_.db_cols;
  lwe::Vector query_vector = lwe::Vector::Zero(params_.db_cols);
  query_vector[col_idx] = 1;

  //RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlace(
  //    query_vector, lwe_pad, log_scaling_factor, lwe_enc_prng.get()));
  
  //RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlaceExplicit(
  //    query_vector, lwe_pad, scaling_factor, lwe_enc_prng.get()));

  RLWE_RETURN_IF_ERROR(s_lwe.EncryptFromPadInPlaceExplicitGivenAs(
      query_vector, As, scaling_factor, lwe_enc_prng.get()));
  

  // Cache the per request state.
  state_ = ClientState{.row_idx = row_idx,
                       .col_idx = col_idx,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};
  state_new_ = ClientStateNew{.row_idx_begin = row_idx,
                        .row_idx_end = row_idx_end,
                       .col_idx = col_idx,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};

  HintlessPirRequest request;
  *request.mutable_ct_query_vector() = SerializeLweCiphertext(query_vector);

  double start, end;
  start = currentDateTime();
  // Step 2. Encrypting the LWE secret using LinPir.
  RLWE_RETURN_IF_ERROR(
      GenerateLinPirRequestInPlace(request, s_lwe.Key()));
  end = currentDateTime();
  std::cout << "[==> TIMER  <==] Client online (LinPir) RLWE-only (rotation keys, etc.) req gen time: " << (end-start) << " ms | " << (end-start)/1000 << " sec" << std::endl;

  return std::make_pair(request, query_vector);
}

absl::StatusOr<std::pair<HintlessPirRequest, lwe::Vector>> Client::GenerateRequestBSGS(int64_t index, lwe::Vector& lwe_sk_vsp_out) {
  if (index < 0 || index >= params_.db_rows * params_.db_cols) {
    return absl::InvalidArgumentError("`index` out of range.");
  }

  // Step 1. Encrypting the selection vector under LWE.
  lwe::Matrix lwe_pad;
  std::unique_ptr<rlwe::SecurePrng> lwe_enc_prng;
  std::string prng_seed_linpir_sk;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(auto pad_prng, rlwe::SingleThreadHkdfPrng::Create(
                                             prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        lwe_pad, lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim,
                                pad_prng.get()));
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadHkdfPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());

  } else {
    RLWE_ASSIGN_OR_RETURN(auto pad_prng, rlwe::SingleThreadChaChaPrng::Create(
                                             prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        lwe_pad, lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim,
                                pad_prng.get()));
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadChaChaPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
  }

  /*
  // Reduce pad modulo the prime (if prime used)
  uint64_t rows = lwe_pad.rows();
  uint64_t cols = lwe_pad.cols();
  for (uint64_t i = 0; i < rows; i++) {
    for (uint64_t j = 0; j < cols; j++) {
      lwe_pad(i, j) = fast_mod_max_64(lwe_pad(i, j));
    }
  }
  */

  RLWE_ASSIGN_OR_RETURN(
      lwe::SymmetricLweKey lwe_secret_key,
      lwe::SymmetricLweKey::Sample(params_.lwe_secret_dim, lwe_enc_prng.get()));

  // Choosing the largest scaling factor that supports our plaintext space
  int log_scaling_factor =
      params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
  lwe::Integer scaling_factor = static_cast<lwe::Integer>(SOL_NTT_PRIME)/(1 << params_.lwe_plaintext_bit_size);
  if (scaling_factor == 0) {
    scaling_factor = (1 << log_scaling_factor);
  }

  // Plaintext is a selection vector for col_idx
  int64_t row_idx = index / params_.db_cols;
  uint64_t row_idx_end = row_idx + params_.db_stack_cells;
  int64_t col_idx = index % params_.db_cols;
  lwe::Vector query_vector = lwe::Vector::Zero(params_.db_cols);
  query_vector[col_idx] = 1;

  //RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlace(
  //    query_vector, lwe_pad, log_scaling_factor, lwe_enc_prng.get()));
  
  RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlaceExplicit(
      query_vector, lwe_pad, scaling_factor, lwe_enc_prng.get()));
  

  // Cache the per request state.
  state_ = ClientState{.row_idx = row_idx,
                       .col_idx = col_idx,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};
  state_new_ = ClientStateNew{.row_idx_begin = row_idx,
                        .row_idx_end = row_idx_end,
                       .col_idx = col_idx,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};

  HintlessPirRequest request;
  *request.mutable_ct_query_vector() = SerializeLweCiphertext(query_vector);

  double start, end;
  start = currentDateTime();
  // Step 2. Encrypting the LWE secret using LinPir.
  RLWE_RETURN_IF_ERROR(
      GenerateLinPirRequestInPlaceBSGS(request, lwe_secret_key.Key()));
  end = currentDateTime();
  std::cout << "[==> TIMER  <==] Client online (LinPir) RLWE-only (BSGS, rotation keys, etc.) req gen time: " << (end-start) << " ms | " << (end-start)/1000 << " sec" << std::endl;

  // Send the LWE SK to VSPIR
  for (int i = 0; i < lwe_secret_key.Key().size(); i++) {
    //std::cout << "LWE SK[" << i << "]: " << lwe_secret_key.Key()(i) << std::endl;
    lwe_sk_vsp_out(i) = lwe_secret_key.Key()(i);
    //std::cout << "LWE SK OUT[" << i << "]: " << key_out(i) << std::endl;
  }
  return std::make_pair(request, query_vector);
}

absl::StatusOr<std::pair<HintlessPirRequest, lwe::Vector>> Client::GenerateRequestBSGSGivenAs(int64_t index, lwe::Vector& As, lwe::SymmetricLweKey& s_lwe) {
  if (index < 0 || index >= params_.db_rows * params_.db_cols) {
    return absl::InvalidArgumentError("`index` out of range.");
  }

  // Step 1. Encrypting the selection vector under LWE.
  std::unique_ptr<rlwe::SecurePrng> lwe_enc_prng;
  std::string prng_seed_linpir_sk;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadHkdfPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());

  } else {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadChaChaPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
  }

  // Choosing the largest scaling factor that supports our plaintext space
  int log_scaling_factor =
      params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
  lwe::Integer scaling_factor = static_cast<lwe::Integer>(SOL_NTT_PRIME)/(1 << params_.lwe_plaintext_bit_size);
  if (scaling_factor == 0) {
    scaling_factor = (1 << log_scaling_factor);
  }

  // Plaintext is a selection vector for col_idx
  int64_t row_idx = index / params_.db_cols;
  uint64_t row_idx_end = row_idx + params_.db_stack_cells;
  int64_t col_idx = index % params_.db_cols;
  lwe::Vector query_vector = lwe::Vector::Zero(params_.db_cols);
  query_vector[col_idx] = 1;

  //RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlace(
  //    query_vector, lwe_pad, log_scaling_factor, lwe_enc_prng.get()));
  
  //RLWE_RETURN_IF_ERROR(s_lwe.EncryptFromPadInPlaceExplicit(
  //    query_vector, lwe_pad, scaling_factor, lwe_enc_prng.get()));

  RLWE_RETURN_IF_ERROR(s_lwe.EncryptFromPadInPlaceExplicitGivenAs(
      query_vector, As, scaling_factor, lwe_enc_prng.get()));
  

  // Cache the per request state.
  state_ = ClientState{.row_idx = row_idx,
                       .col_idx = col_idx,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};
  state_new_ = ClientStateNew{.row_idx_begin = row_idx,
                        .row_idx_end = row_idx_end,
                       .col_idx = col_idx,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};

  HintlessPirRequest request;
  *request.mutable_ct_query_vector() = SerializeLweCiphertext(query_vector);

  double start, end;
  start = currentDateTime();
  // Step 2. Encrypting the LWE secret using LinPir.
  RLWE_RETURN_IF_ERROR(
      GenerateLinPirRequestInPlaceBSGS(request, s_lwe.Key()));
  end = currentDateTime();
  std::cout << BLUE << "[==> TIMER  <==] Client online (LinPir) RLWE-only (BSGS, rotation keys, etc.) req gen time: " << (end-start) << " ms | " << (end-start)/1000 << " sec" << END << std::endl;

  return std::make_pair(request, query_vector);
}

absl::StatusOr<HintlessPirRequest> Client::PrepareLinPirGivenS_BSGS(lwe::SymmetricLweKey& s_lwe) {

  std::string prng_seed_linpir_sk;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());

  } else {
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
  }
  

  // Cache the per request state.
  state_ = ClientState{.row_idx = 0,
                       .col_idx = 0,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};
  state_new_ = ClientStateNew{.row_idx_begin = 0,
                        .row_idx_end = 0,
                       .col_idx = 0,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};

  HintlessPirRequest prepare_request;

  RLWE_RETURN_IF_ERROR(
      GenerateLinPirRequestInPlaceBSGS(prepare_request, s_lwe.Key()));

  return prepare_request;
}

absl::StatusOr<HintlessPirRequest> Client::PrepareLinPirGivenS(lwe::SymmetricLweKey& s_lwe) {

  std::string prng_seed_linpir_sk;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());

  } else {
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
  }
  

  // Cache the per request state.
  state_ = ClientState{.row_idx = 0,
                       .col_idx = 0,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};
  state_new_ = ClientStateNew{.row_idx_begin = 0,
                        .row_idx_end = 0,
                       .col_idx = 0,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};

  HintlessPirRequest prepare_request;

  RLWE_RETURN_IF_ERROR(
      GenerateLinPirRequestInPlace(prepare_request, s_lwe.Key()));

  return prepare_request;
}

absl::StatusOr<std::pair<HintlessPirRequest, lwe::Vector>> Client::GenerateRequestGivenAsSkipLinPir(int64_t index, lwe::Vector& As, lwe::SymmetricLweKey& s_lwe) {
  if (index < 0 || index >= params_.db_rows * params_.db_cols) {
    return absl::InvalidArgumentError("`index` out of range.");
  }

  // Step 1. Encrypting the selection vector under LWE.
  std::unique_ptr<rlwe::SecurePrng> lwe_enc_prng;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadHkdfPrng::Create(prng_seed_enc));

  } else {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadChaChaPrng::Create(prng_seed_enc));
  }

  // Choosing the largest scaling factor that supports our plaintext space
  int log_scaling_factor =
      params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
  lwe::Integer scaling_factor = static_cast<lwe::Integer>(SOL_NTT_PRIME)/(1 << params_.lwe_plaintext_bit_size);
  if (scaling_factor == 0) {
    scaling_factor = (1 << log_scaling_factor);
  }

  // Plaintext is a selection vector for col_idx
  int64_t row_idx = index / params_.db_cols;
  uint64_t row_idx_end = row_idx + params_.db_stack_cells;
  int64_t col_idx = index % params_.db_cols;
  lwe::Vector query_vector = lwe::Vector::Zero(params_.db_cols);
  query_vector[col_idx] = 1;

  //RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlace(
  //    query_vector, lwe_pad, log_scaling_factor, lwe_enc_prng.get()));
  
  //RLWE_RETURN_IF_ERROR(s_lwe.EncryptFromPadInPlaceExplicit(
  //    query_vector, lwe_pad, scaling_factor, lwe_enc_prng.get()));

  RLWE_RETURN_IF_ERROR(s_lwe.EncryptFromPadInPlaceExplicitGivenAs(
      query_vector, As, scaling_factor, lwe_enc_prng.get()));
  

  // Cache the per request state.
  state_.row_idx = row_idx;
  state_.col_idx = col_idx;
                       
  state_new_.row_idx_begin = row_idx;
  state_new_.row_idx_end = row_idx_end;
  state_new_.col_idx = col_idx;

  HintlessPirRequest request;
  *request.mutable_ct_query_vector() = SerializeLweCiphertext(query_vector);

  return std::make_pair(request, query_vector);
}

absl::StatusOr<std::tuple<VMatrix, lwe::Vector, lwe::SymmetricLweKey>> Client::Compute_A_times_s() {
  double start, end;
  start = currentDateTime();
  // Step 1. Encrypting the selection vector under LWE.
  lwe::Matrix lwe_pad;
  std::unique_ptr<rlwe::SecurePrng> lwe_enc_prng;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(auto pad_prng, rlwe::SingleThreadHkdfPrng::Create(
                                             prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        lwe_pad, lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim,
                                pad_prng.get()));
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadHkdfPrng::Create(prng_seed_enc));

  } else {
    RLWE_ASSIGN_OR_RETURN(auto pad_prng, rlwe::SingleThreadChaChaPrng::Create(
                                             prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        lwe_pad, lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim,
                                pad_prng.get()));
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadChaChaPrng::Create(prng_seed_enc));
  }
  end = currentDateTime();
  //std::cout << BLUE << "[==> TIMER  <==] Client expanding A from seed time: " << (end-start) << " ms | " << (end-start)/1000 << " sec" << END << std::endl;

  start = currentDateTime();

 	VMatrix A_vsp(params_.db_cols, params_.lwe_secret_dim);
  lwe::Vector sk_vsp = lwe::Vector::Zero(params_.lwe_secret_dim);
  lwe::Vector As = lwe::Vector::Zero(params_.db_cols);

  convertToVeriSimplePIRMatrix(A_vsp, lwe_pad);
  //std::cout << YELLOW << "[ONLINE STATE]: Size of of A[0] (first row) " << (A_vsp.cols * sizeof(Elem))/(1ULL << 10) << " KB" << END << std::endl;

  RLWE_ASSIGN_OR_RETURN(
      lwe::SymmetricLweKey lwe_secret_key,
      lwe::SymmetricLweKey::Sample(params_.lwe_secret_dim, lwe_enc_prng.get()));

  for (int i = 0; i < lwe_secret_key.Key().size(); i++) {
    sk_vsp(i) = lwe_secret_key.Key()(i);
  }
  VMatrix sk_vsp__ = convertToVeriSimplePIRColumnMatrix(sk_vsp);
  auto As_vsp = matMulVec(A_vsp, sk_vsp__);
  convertToVectorFromVeriSimplePIRColumnMatrix(As, As_vsp);
  
  end = currentDateTime();
  //std::cout << BLUE << "[==> TIMER  <==] Client A*s time: " << (end-start) << " ms | " << (end-start)/1000 << " sec" << END << std::endl;

  return std::make_tuple(As_vsp, As, lwe_secret_key);
}

absl::StatusOr<HintlessPirRequest> Client::GenerateRequestWithoutQuery(int64_t index, lwe::Vector& lwe_sk_vsp_out, lwe::Vector& query_vsp_out) {
  if (index < 0 || index >= params_.db_rows * params_.db_cols) {
    return absl::InvalidArgumentError("`index` out of range.");
  }

  // Step 1. Encrypting the selection vector under LWE.
  lwe::Matrix lwe_pad;
  std::unique_ptr<rlwe::SecurePrng> lwe_enc_prng;
  std::string prng_seed_linpir_sk;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(auto pad_prng, rlwe::SingleThreadHkdfPrng::Create(
                                             prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        lwe_pad, lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim,
                                pad_prng.get()));
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadHkdfPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());

  } else {
    RLWE_ASSIGN_OR_RETURN(auto pad_prng, rlwe::SingleThreadChaChaPrng::Create(
                                             prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        lwe_pad, lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim,
                                pad_prng.get()));
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed_enc,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(lwe_enc_prng,
                          rlwe::SingleThreadChaChaPrng::Create(prng_seed_enc));
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_sk,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
  }

  // Reduce pad modulo the prime (if prime used)
  uint64_t rows = lwe_pad.rows();
  uint64_t cols = lwe_pad.cols();
  for (uint64_t i = 0; i < rows; i++) {
    for (uint64_t j = 0; j < cols; j++) {
      lwe_pad(i, j) = fast_mod_max_64(lwe_pad(i, j));
    }
  }

  RLWE_ASSIGN_OR_RETURN(
      lwe::SymmetricLweKey lwe_secret_key,
      lwe::SymmetricLweKey::Sample(params_.lwe_secret_dim, lwe_enc_prng.get()));

  // Choosing the largest scaling factor that supports our plaintext space
  int log_scaling_factor =
      params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;

  // Plaintext is a selection vector for col_idx
  int64_t row_idx = index / params_.db_cols;
  int64_t col_idx = index % params_.db_cols;
  //lwe::Vector query_vector = lwe::Vector::Zero(params_.db_cols);
  //query_vector[col_idx] = 1;
  query_vsp_out[col_idx] = 1;

  //RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlace(
  //    query_vector, lwe_pad, log_scaling_factor, lwe_enc_prng.get()));
  RLWE_RETURN_IF_ERROR(lwe_secret_key.EncryptFromPadInPlace(
      query_vsp_out, lwe_pad, log_scaling_factor, lwe_enc_prng.get()));

  // Cache the per request state.
  state_ = ClientState{.row_idx = row_idx,
                       .col_idx = col_idx,
                       .prng_seed_linpir_sk = std::move(prng_seed_linpir_sk)};

  HintlessPirRequest request;
  // Query will be sent by VSPIR, so we ignore sending from here
  //*request.mutable_ct_query_vector() = SerializeLweCiphertext(query_vector);

  // Step 2. Encrypting the LWE secret using LinPir.
  RLWE_RETURN_IF_ERROR(
      GenerateLinPirRequestInPlace(request, lwe_secret_key.Key()));

  // Send the LWE SK to VSPIR
  for (int i = 0; i < lwe_secret_key.Key().size(); i++) {
    //std::cout << "LWE SK[" << i << "]: " << lwe_secret_key.Key()(i) << std::endl;
    lwe_sk_vsp_out(i) = lwe_secret_key.Key()(i);
    //std::cout << "LWE SK OUT[" << i << "]: " << key_out(i) << std::endl;
  }
  
  return request;
}

absl::Status Client::GenerateLinPirRequestInPlace(
    HintlessPirRequest& request, const lwe::Vector& lwe_secret) const {
  if (linpir_clients_.empty()) {
    return absl::InvalidArgumentError("No LinPir client available.");
  }

  // Encode the LWE secret vector using LinPir plaintext moduli, and also
  // generate a GaloisKey which is shared by all LinPir requests.
  RlweInteger lwe_modulus = RlweInteger{1} << params_.lwe_modulus_bit_size;
  if (SOL_NTT_PRIME > 0) {
    lwe_modulus = SOL_NTT_PRIME;
  }
  //std::cout << "THIS THIS " << lwe_modulus << std::endl;
  for (int k = 0; k < linpir_clients_.size(); ++k) {
    RlweInteger plaintext_modulus = rlwe_contexts_[k]->PlaintextModulus();
    std::vector<RlweInteger> lwe_secret_mod_t =
        EncodeLweVector(lwe_secret, lwe_modulus, plaintext_modulus);
    RLWE_ASSIGN_OR_RETURN(
        auto ct, linpir_clients_[k]->EncryptQuery(lwe_secret_mod_t,
                                                  state_.prng_seed_linpir_sk));
    RLWE_ASSIGN_OR_RETURN(auto ct_b, ct.Component(0));
    RLWE_ASSIGN_OR_RETURN(*request.add_linpir_ct_bs(),
                          ct_b.Serialize(rlwe_moduli_));
  }
  RLWE_ASSIGN_OR_RETURN(auto gk, linpir_clients_[0]->GenerateGaloisKey(
                                     state_.prng_seed_linpir_sk));
  for (auto const& gk_b : gk.GetKeyB()) {
    RLWE_ASSIGN_OR_RETURN(*request.add_linpir_gk_bs(),
                          gk_b.Serialize(rlwe_moduli_));
  }
  return absl::OkStatus();
}

absl::Status Client::GenerateLinPirRequestInPlaceBSGS(
    HintlessPirRequest& request, const lwe::Vector& lwe_secret) const {
  if (linpir_clients_.empty()) {
    return absl::InvalidArgumentError("No LinPir client available.");
  }

  // Encode the LWE secret vector using LinPir plaintext moduli, and also
  // generate a GaloisKey which is shared by all LinPir requests.
  RlweInteger lwe_modulus = RlweInteger{1} << params_.lwe_modulus_bit_size;
  if (SOL_NTT_PRIME > 0) {
    lwe_modulus = SOL_NTT_PRIME;
  }
  //std::cout << "THIS THIS " << lwe_modulus << std::endl;
  for (int k = 0; k < linpir_clients_.size(); ++k) {
    RlweInteger plaintext_modulus = rlwe_contexts_[k]->PlaintextModulus();
    std::vector<RlweInteger> lwe_secret_mod_t =
        EncodeLweVector(lwe_secret, lwe_modulus, plaintext_modulus);
    RLWE_ASSIGN_OR_RETURN(
        auto ct, linpir_clients_[k]->EncryptQuery(lwe_secret_mod_t,
                                                  state_.prng_seed_linpir_sk));
    RLWE_ASSIGN_OR_RETURN(auto ct_b, ct.Component(0));
    RLWE_ASSIGN_OR_RETURN(*request.add_linpir_ct_bs(),
                          ct_b.Serialize(rlwe_moduli_));
  }
  RLWE_ASSIGN_OR_RETURN(auto gk, linpir_clients_[0]->GenerateGaloisKeyBSGS(
                                     state_.prng_seed_linpir_sk));
  //auto gk = linpir_clients_[0]->GenerateGaloisKeyBSGS(state_.prng_seed_linpir_sk)
  auto gk_baby = gk.first;
  auto gk_giant = gk.second;
  for (auto const& gk_b : gk_baby.GetKeyB()) {
    RLWE_ASSIGN_OR_RETURN(*request.add_linpir_gk_bs_baby(),
                          gk_b.Serialize(rlwe_moduli_));
  }
  for (auto const& gk_b : gk_giant.GetKeyB()) {
    RLWE_ASSIGN_OR_RETURN(*request.add_linpir_gk_bs_giant(),
                          gk_b.Serialize(rlwe_moduli_));
  }
  return absl::OkStatus();
} 

std::vector<Client::RlweInteger> Client::EncodeLweVector(
    const lwe::Vector& lwe_vector, RlweInteger lwe_modulus,
    RlweInteger encode_modulus) {
  RlweInteger lwe_modulus_half = lwe_modulus >> 1;
  std::vector<RlweInteger> lwe_vector_mod_t(lwe_vector.size(), 0);
  for (int i = 0; i < lwe_vector.size(); ++i) {
    RlweInteger x = lwe_vector[i];
    lwe_vector_mod_t[i] =
        ConvertModulus(x, lwe_modulus, encode_modulus, lwe_modulus_half);
  }
  return lwe_vector_mod_t;
}

absl::StatusOr<std::string> Client::RecoverRecord(
    const HintlessPirResponse& response) {
  int num_shards =
      DivAndRoundUp(params_.db_record_bit_size, params_.lwe_plaintext_bit_size);
  if (response.ct_records_size() != num_shards) {
    return absl::InvalidArgumentError("`response` has incorrect size.");
  }

  // Recover decryption_parts = Hint * LWE secret = Database * A * LWE secret.
  RLWE_ASSIGN_OR_RETURN(std::vector<lwe::Vector> decryption_parts,
                        RecoverLweDecryptionParts(response));

  // Decrypt the LWE ciphertexts in response.
  std::vector<lwe::Integer> values;
  values.reserve(response.ct_records_size());
  for (int i = 0; i < response.ct_records_size(); ++i) {
    
    std::vector<lwe::Integer> ct_records =
        DeserializeLweCiphertext(response.ct_records(i));
    if (ct_records.size() != params_.db_rows) {
      return absl::InvalidArgumentError(absl::StrCat(
          "The server response has incorrect dimension; got ",
          ct_records.size(), " but expecting ", params_.db_rows, "."));
    }

    // Remove hint * s from the server response, which gives us \Delta * m + e.
    lwe::Vector noisy_plaintext{{ct_records[state_.row_idx]}};
    
    noisy_plaintext[0] -= decryption_parts[i][state_.row_idx];

    // Remove the error e.
    int log_scaling_factor =
        params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
    RLWE_RETURN_IF_ERROR(
        lwe::RemoveErrorInPlace(noisy_plaintext, log_scaling_factor));

    // Extracting the coefficient from the 1 x 1 matrix noisy_plaintext.
    values.push_back(noisy_plaintext.eval()(0));
  }

  return ReconstructRecord(values, params_);
}

absl::StatusOr<std::tuple<std::string, std::vector<lwe::Vector>, std::vector<std::vector<lwe::Integer>>>> Client::RecoverRecordAndGiveHsAndV(
    const HintlessPirResponse& response) {
  int num_shards =
      DivAndRoundUp(params_.db_record_bit_size, params_.lwe_plaintext_bit_size);

  auto plaintext_moduli = crt_context_.MainPrimeModuli();
  int response_size = response.linpir_responses_size() / plaintext_moduli.size();
  //std::cout << "response size " << response_size << std::endl;
  //std::cout << "num shards " << num_shards << std::endl;
  
  if (response.ct_records_size() != num_shards) {
    return absl::InvalidArgumentError("`response` has incorrect size.");
  }

  //std::cout << "ct records size " << response.ct_records_size() << std::endl;
  //std::cout << "linpir responses size " << response.linpir_responses_size() << std::endl;
  //std::cout << "plaintext moduli size " << plaintext_moduli.size() << std::endl;

  // Recover decryption_parts = Hint * LWE secret = Database * A * LWE secret.
  RLWE_ASSIGN_OR_RETURN(std::vector<lwe::Vector> decryption_parts,
                        RecoverLweDecryptionParts(response));

  //std::cout << "H*s recovered from RLWE" << std::endl;

  // Decrypt the LWE ciphertexts in response.
  std::vector<lwe::Integer> values;
  std::vector<std::vector<lwe::Integer>> pir_response_ct_vec(num_shards);
  values.reserve(response.ct_records_size() * params_.db_stack_cells);
  //pir_response_ct_vec.reserve(num_shards);
  
  //assert(response.ct_records_size() == 1 && "Temporarily, only one decryption part is supported");
  //std::cout << "A" << std::endl;
  for (int i = 0; i < response.ct_records_size(); ++i) {
    //pir_response_ct_vec[i].reserve(params_.db_rows);
    std::vector<lwe::Integer> ct_records =
        DeserializeLweCiphertext(response.ct_records(i));
    if (ct_records.size() != params_.db_rows) {
      return absl::InvalidArgumentError(absl::StrCat(
          "The server response has incorrect dimension; got ",
          ct_records.size(), " but expecting ", params_.db_rows, "."));
    }
    //std::cout << "B" << std::endl;
    // initialize pir_response_ct_vector[i]
    pir_response_ct_vec[i] = std::vector<lwe::Integer>(ct_records.size());

    for (int j = 0; j < ct_records.size(); j++) {
      pir_response_ct_vec[i][j] = ct_records[j];
    }

    // Remove hint * s from the server response, which gives us \Delta * m + e.
    // Vector of size 1

    if (params_.db_stack_cells == 1){
      lwe::Vector noisy_plaintext{{ct_records[state_.row_idx]}};

      // We only care about decrypting a specific entry, so we can ignore the rest
      //noisy_plaintext[0] -= decryption_parts[i][state_.row_idx];
      
      Elem64 tmp = ((Elem64)noisy_plaintext[0]) + ((Elem64)SOL_NTT_PRIME) - ((Elem64)decryption_parts[i][state_.row_idx]);
      tmp = fast_mod_max_64(tmp);
      noisy_plaintext[0] = static_cast<lwe::Integer>(tmp);

      // Remove the error e.
      int log_scaling_factor =
          params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
      lwe::Integer scaling_factor = static_cast<lwe::Integer>(SOL_NTT_PRIME)/(1 << params_.lwe_plaintext_bit_size);
      if (scaling_factor == 0) {
        scaling_factor = (1 << log_scaling_factor);
      }
      // TODO: Instead of log_scaling_factor, have full scaling factor. Rewrite RemoveErrorInPlace.
      //RLWE_RETURN_IF_ERROR(
      //    lwe::RemoveErrorInPlace(noisy_plaintext, log_scaling_factor));
      RLWE_RETURN_IF_ERROR(
          lwe::RemoveErrorInPlaceCustomMod(noisy_plaintext, scaling_factor, (1 << params_.lwe_plaintext_bit_size)));

      // Extracting the coefficient from the 1 x 1 matrix noisy_plaintext.
      values.push_back(noisy_plaintext.eval()(0)); // eval means lazy evaluation
    } 
    else { // When stacking is implemented
      //lwe::Vector noisy_plaintext;
      //noisy_plaintext.reserve(params_.db_stack_cells);
      for (int j = 0; j < params_.db_stack_cells; j++) {
        //std::cout << i << " " << j << std::endl;
        lwe::Vector noisy_plaintext{{ct_records[state_new_.row_idx_begin + j]}};

        // We only care about decrypting a specific entry, so we can ignore the rest
        //noisy_plaintext[0] -= decryption_parts[i][state_.row_idx];
        
        Elem64 tmp = ((Elem64)noisy_plaintext[0]) + ((Elem64)SOL_NTT_PRIME) - ((Elem64)decryption_parts[i][state_new_.row_idx_begin + j]);
        tmp = fast_mod_max_64(tmp);
        noisy_plaintext[0] = static_cast<lwe::Integer>(tmp);

        // Remove the error e.
        int log_scaling_factor =
            params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
        lwe::Integer scaling_factor = static_cast<lwe::Integer>(SOL_NTT_PRIME)/(1 << params_.lwe_plaintext_bit_size);
        if (scaling_factor == 0) {
          scaling_factor = (1 << log_scaling_factor);
        }
        // TODO: Instead of log_scaling_factor, have full scaling factor. Rewrite RemoveErrorInPlace.
        //RLWE_RETURN_IF_ERROR(
        //    lwe::RemoveErrorInPlace(noisy_plaintext, log_scaling_factor));
        RLWE_RETURN_IF_ERROR(
            lwe::RemoveErrorInPlaceCustomMod(noisy_plaintext, scaling_factor, (1 << params_.lwe_plaintext_bit_size)));

        // Extracting the coefficient from the 1 x 1 matrix noisy_plaintext.
        values.push_back(noisy_plaintext.eval()(0)); // eval means lazy evaluation
      }
    }
  }
  //assert(decryption_parts.size() == 1 && "Temporarily, only one decryption part is supported");
  std::cout << "Reconstructing record" << std::endl;
  if(params_.db_stack_cells > 1) {
    return std::make_tuple(ReconstructRecordWithStacking(values, params_), decryption_parts, pir_response_ct_vec);
  }
  return std::make_tuple(ReconstructRecord(values, params_), decryption_parts, pir_response_ct_vec);
}

absl::StatusOr<std::vector<lwe::Vector>> Client::RecoverHsPreparePhase(
    const HintlessPirResponse& response) {
  // Recover decryption_parts = Hint * LWE secret = Database * A * LWE secret.
  RLWE_ASSIGN_OR_RETURN(std::vector<lwe::Vector> decryption_parts,
                        RecoverLweDecryptionParts(response));

  std::cout << "H*s recovered from RLWE" << std::endl;
  return decryption_parts;
}

absl::StatusOr<std::tuple<std::string, std::vector<std::vector<lwe::Integer>>>> Client::RecoverRecordAndGiveV(
    const HintlessPirResponse& response, const std::vector<lwe::Vector>& decryption_parts) {
  int num_shards =
      DivAndRoundUp(params_.db_record_bit_size, params_.lwe_plaintext_bit_size);

  auto plaintext_moduli = crt_context_.MainPrimeModuli();
  int response_size = response.linpir_responses_size() / plaintext_moduli.size();
  //std::cout << "response size " << response_size << std::endl;
  //std::cout << "num shards " << num_shards << std::endl;
  
  if (response.ct_records_size() != num_shards) {
    return absl::InvalidArgumentError("`response` has incorrect size.");
  }

  //std::cout << "ct records size " << response.ct_records_size() << std::endl;
  //std::cout << "linpir responses size " << response.linpir_responses_size() << std::endl;
  //std::cout << "plaintext moduli size " << plaintext_moduli.size() << std::endl;

  // Decrypt the LWE ciphertexts in response.
  std::vector<lwe::Integer> values;
  std::vector<std::vector<lwe::Integer>> pir_response_ct_vec(num_shards);
  values.reserve(response.ct_records_size() * params_.db_stack_cells);
  //pir_response_ct_vec.reserve(num_shards);
  
  //assert(response.ct_records_size() == 1 && "Temporarily, only one decryption part is supported");
  //std::cout << "A" << std::endl;
  for (int i = 0; i < response.ct_records_size(); ++i) {
    //pir_response_ct_vec[i].reserve(params_.db_rows);
    std::vector<lwe::Integer> ct_records =
        DeserializeLweCiphertext(response.ct_records(i));
    if (ct_records.size() != params_.db_rows) {
      return absl::InvalidArgumentError(absl::StrCat(
          "The server response has incorrect dimension; got ",
          ct_records.size(), " but expecting ", params_.db_rows, "."));
    }
    //std::cout << "B" << std::endl;
    // initialize pir_response_ct_vector[i]
    pir_response_ct_vec[i] = std::vector<lwe::Integer>(ct_records.size());

    for (int j = 0; j < ct_records.size(); j++) {
      pir_response_ct_vec[i][j] = ct_records[j];
    }

    // Remove hint * s from the server response, which gives us \Delta * m + e.
    // Vector of size 1

    if (params_.db_stack_cells == 1){
      lwe::Vector noisy_plaintext{{ct_records[state_.row_idx]}};

      // We only care about decrypting a specific entry, so we can ignore the rest
      //noisy_plaintext[0] -= decryption_parts[i][state_.row_idx];
      
      Elem64 tmp = ((Elem64)noisy_plaintext[0]) + ((Elem64)SOL_NTT_PRIME) - ((Elem64)decryption_parts[i][state_.row_idx]);
      tmp = fast_mod_max_64(tmp);
      noisy_plaintext[0] = static_cast<lwe::Integer>(tmp);

      // Remove the error e.
      int log_scaling_factor =
          params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
      lwe::Integer scaling_factor = static_cast<lwe::Integer>(SOL_NTT_PRIME)/(1 << params_.lwe_plaintext_bit_size);
      if (scaling_factor == 0) {
        scaling_factor = (1 << log_scaling_factor);
      }
      // TODO: Instead of log_scaling_factor, have full scaling factor. Rewrite RemoveErrorInPlace.
      //RLWE_RETURN_IF_ERROR(
      //    lwe::RemoveErrorInPlace(noisy_plaintext, log_scaling_factor));
      RLWE_RETURN_IF_ERROR(
          lwe::RemoveErrorInPlaceCustomMod(noisy_plaintext, scaling_factor, (1 << params_.lwe_plaintext_bit_size)));

      // Extracting the coefficient from the 1 x 1 matrix noisy_plaintext.
      values.push_back(noisy_plaintext.eval()(0)); // eval means lazy evaluation
    } 
    else { // When stacking is implemented
      //lwe::Vector noisy_plaintext;
      //noisy_plaintext.reserve(params_.db_stack_cells);
      for (int j = 0; j < params_.db_stack_cells; j++) {
        //std::cout << i << " " << j << std::endl;
        lwe::Vector noisy_plaintext{{ct_records[state_new_.row_idx_begin + j]}};

        // We only care about decrypting a specific entry, so we can ignore the rest
        //noisy_plaintext[0] -= decryption_parts[i][state_.row_idx];
        
        Elem64 tmp = ((Elem64)noisy_plaintext[0]) + ((Elem64)SOL_NTT_PRIME) - ((Elem64)decryption_parts[i][state_new_.row_idx_begin + j]);
        tmp = fast_mod_max_64(tmp);
        noisy_plaintext[0] = static_cast<lwe::Integer>(tmp);

        // Remove the error e.
        int log_scaling_factor =
            params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
        lwe::Integer scaling_factor = static_cast<lwe::Integer>(SOL_NTT_PRIME)/(1 << params_.lwe_plaintext_bit_size);
        if (scaling_factor == 0) {
          scaling_factor = (1 << log_scaling_factor);
        }
        // TODO: Instead of log_scaling_factor, have full scaling factor. Rewrite RemoveErrorInPlace.
        //RLWE_RETURN_IF_ERROR(
        //    lwe::RemoveErrorInPlace(noisy_plaintext, log_scaling_factor));
        RLWE_RETURN_IF_ERROR(
            lwe::RemoveErrorInPlaceCustomMod(noisy_plaintext, scaling_factor, (1 << params_.lwe_plaintext_bit_size)));

        // Extracting the coefficient from the 1 x 1 matrix noisy_plaintext.
        values.push_back(noisy_plaintext.eval()(0)); // eval means lazy evaluation
      }
    }
  }
  //assert(decryption_parts.size() == 1 && "Temporarily, only one decryption part is supported");
  std::cout << "Reconstructing record" << std::endl;
  if(params_.db_stack_cells > 1) {
    return std::make_tuple(ReconstructRecordWithStacking(values, params_), pir_response_ct_vec);
  }
  return std::make_tuple(ReconstructRecord(values, params_), pir_response_ct_vec);
}

absl::StatusOr<std::tuple<std::string, lwe::Vector, std::vector<lwe::Integer>>> Client::RecoverRecordAndGiveHsAndVNaive(
    const HintlessPirResponse& response, std::vector<lwe::Vector>& decryption_parts) {
  int num_shards =
      DivAndRoundUp(params_.db_record_bit_size, params_.lwe_plaintext_bit_size);

  auto plaintext_moduli = crt_context_.MainPrimeModuli();
  int response_size = response.linpir_responses_size() / plaintext_moduli.size();
  //std::cout << "response size " << response_size << std::endl;
  //std::cout << "num shards " << num_shards << std::endl;
  
  if (response.ct_records_size() != num_shards) {
    return absl::InvalidArgumentError("`response` has incorrect size.");
  }

  //std::cout << "ct records size " << response.ct_records_size() << std::endl;
  //std::cout << "linpir responses size " << response.linpir_responses_size() << std::endl;
  //std::cout << "plaintext moduli size " << plaintext_moduli.size() << std::endl;

  // Recover decryption_parts = Hint * LWE secret = Database * A * LWE secret.
  //RLWE_ASSIGN_OR_RETURN(std::vector<lwe::Vector> decryption_parts,
  //                      RecoverLweDecryptionParts(response));

  //std::cout << "H*s recovered from RLWE" << std::endl;

  // Decrypt the LWE ciphertexts in response.
  std::vector<lwe::Integer> values, pir_response_ct;
  values.reserve(response.ct_records_size());
  pir_response_ct.reserve(params_.db_rows);
  assert(response.ct_records_size() == 1 && "Temporarily, only one decryption part is supported");
  for (int i = 0; i < response.ct_records_size(); ++i) {
    
    std::vector<lwe::Integer> ct_records =
        DeserializeLweCiphertext(response.ct_records(i));
    if (ct_records.size() != params_.db_rows) {
      return absl::InvalidArgumentError(absl::StrCat(
          "The server response has incorrect dimension; got ",
          ct_records.size(), " but expecting ", params_.db_rows, "."));
    }

    for (int j = 0; j < ct_records.size(); j++) {
      pir_response_ct.push_back(ct_records[j]);
    }
    

    // Remove hint * s from the server response, which gives us \Delta * m + e.
    // Vector of size 1
    lwe::Vector noisy_plaintext{{ct_records[state_.row_idx]}};
    
    // TODO: Make this prime modulo friendly so that result is modulo prime
    // We only care about decrypting a specific entry, so we can ignore the rest
    //noisy_plaintext[0] -= decryption_parts[i][state_.row_idx];
    
    
    Elem64 tmp = ((Elem64)noisy_plaintext[0]) + ((Elem64)SOL_NTT_PRIME) - ((Elem64)decryption_parts[i][state_.row_idx]);
    tmp = fast_mod_max_64(tmp);
    noisy_plaintext[0] = static_cast<lwe::Integer>(tmp);
    

    // Remove the error e.
    int log_scaling_factor =
        params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
    lwe::Integer scaling_factor = static_cast<lwe::Integer>(SOL_NTT_PRIME)/(1 << params_.lwe_plaintext_bit_size);
    if (scaling_factor == 0) {
      scaling_factor = (1 << log_scaling_factor);
    }
    // TODO: Instead of log_scaling_factor, have full scaling factor. Rewrite RemoveErrorInPlace.
    //RLWE_RETURN_IF_ERROR(
    //    lwe::RemoveErrorInPlace(noisy_plaintext, log_scaling_factor));
    
    RLWE_RETURN_IF_ERROR(
        lwe::RemoveErrorInPlaceCustomMod(noisy_plaintext, scaling_factor, (1 << params_.lwe_plaintext_bit_size)));
    

    // Extracting the coefficient from the 1 x 1 matrix noisy_plaintext.
    values.push_back(noisy_plaintext.eval()(0));
  }
  assert(decryption_parts.size() == 1 && "Temporarily, only one decryption part is supported");
  return std::make_tuple(ReconstructRecord(values, params_), decryption_parts[0], pir_response_ct);
}

absl::StatusOr<std::pair<std::string, lwe::Vector>> Client::RecoverRecordGivenSimplePIRResp(
    const HintlessPirResponse& response, std::vector<std::vector<lwe::Integer>> res_in) {
  int num_shards =
      DivAndRoundUp(params_.db_record_bit_size, params_.lwe_plaintext_bit_size);

  auto plaintext_moduli = crt_context_.MainPrimeModuli();
  int response_size = response.linpir_responses_size() / plaintext_moduli.size();
  //std::cout << "response size " << response_size << std::endl;
  //std::cout << "num shards " << num_shards << std::endl;
  //if (response.ct_records_size() != num_shards) {
  if (response_size != num_shards) {
    return absl::InvalidArgumentError("`response` has incorrect size.");
  }
  //std::cout << "ct records size " << response.ct_records_size() << std::endl;
  //std::cout << "linpir responses size " << response.linpir_responses_size() << std::endl;
  //std::cout << "plaintext moduli size " << plaintext_moduli.size() << std::endl;
  // Recover decryption_parts = Hint * LWE secret = Database * A * LWE secret.
  RLWE_ASSIGN_OR_RETURN(std::vector<lwe::Vector> decryption_parts,
                        RecoverLweDecryptionParts(response));

  // Decrypt the LWE ciphertexts in response.
  std::vector<lwe::Integer> values;
  //values.reserve(response.ct_records_size());
  //assert(res_in.size() == response.ct_records_size());
  //for (int i = 0; i < response.ct_records_size(); ++i) {
  values.reserve(response_size);
  assert(res_in.size() == response_size);
  for (int i = 0; i < response_size; ++i) {  
    
    //std::vector<lwe::Integer> ct_records =
    //    DeserializeLweCiphertext(response.ct_records(i));
    std::vector<lwe::Integer> ct_records = res_in[i];
    if (ct_records.size() != params_.db_rows) {
      return absl::InvalidArgumentError(absl::StrCat(
          "The server response has incorrect dimension; got ",
          ct_records.size(), " but expecting ", params_.db_rows, "."));
    }

    // Remove hint * s from the server response, which gives us \Delta * m + e.
    lwe::Vector noisy_plaintext{{ct_records[state_.row_idx]}};
    
    noisy_plaintext[0] -= decryption_parts[i][state_.row_idx];
    //std::cout << "Row index " << state_.row_idx << " Col index " << state_.col_idx << std::endl; 

    // Remove the error e.
    int log_scaling_factor =
        params_.lwe_modulus_bit_size - params_.lwe_plaintext_bit_size;
    RLWE_RETURN_IF_ERROR(
        lwe::RemoveErrorInPlace(noisy_plaintext, log_scaling_factor));

    // Extracting the coefficient from the 1 x 1 matrix noisy_plaintext.
    values.push_back(noisy_plaintext.eval()(0));
  }

  assert(decryption_parts.size() == 1 && "Temporarily, only one decryption part is supported");
  return std::make_pair(ReconstructRecord(values, params_), decryption_parts[0]);
}

absl::StatusOr<std::vector<lwe::Vector>> Client::RecoverLweDecryptionParts(
    const HintlessPirResponse& response) const {
  using BigInteger = rlwe::uint256;

  auto plaintext_moduli = crt_context_.MainPrimeModuli();
  int num_linpir_plaintext_moduli = plaintext_moduli.size();
  if (response.linpir_responses_size() != num_linpir_plaintext_moduli) {
    return absl::InvalidArgumentError(
        "`response` contains unexpected number of LinPir responses.");
  }

  int num_shards =
      DivAndRoundUp(params_.db_record_bit_size, params_.lwe_plaintext_bit_size);
  if (num_shards != response.linpir_responses(0).ct_inner_products_size()) {
    return absl::InvalidArgumentError(
        "`response` contains an expected number of shards.");
  }

  // CRT decomposed values of Hint * LWE secrets, where the first dimension
  // is per database shard, then per CRT modulus, and then per block in the
  // shard.
  std::vector<std::vector<std::vector<RlweModularInt>>> hint_crt_values(
      num_shards);
  for (auto& h : hint_crt_values) {
    h.resize(num_linpir_plaintext_moduli);
  }
  for (int k = 0; k < num_linpir_plaintext_moduli; ++k) {
    RLWE_ASSIGN_OR_RETURN(
        auto hint_values_mod_tk,
        linpir_clients_[k]->Recover(response.linpir_responses(k)));
    auto mod_params_tk = plaintext_moduli[k]->ModParams();
    for (int i = 0; i < num_shards; ++i) {
      hint_crt_values[i][k].reserve(hint_values_mod_tk[i].size());
      for (auto const& hint_value : hint_values_mod_tk[i]) {
        RLWE_ASSIGN_OR_RETURN(auto hint_mod_tk, RlweModularInt::ImportInt(
                                                    hint_value, mod_params_tk));
        hint_crt_values[i][k].push_back(std::move(hint_mod_tk));
      }
    }
  }
  //std::cout << "Hint CRT values size " << hint_crt_values.size() << " outer CRT slots [2 most likely] " << hint_crt_values[0].size() << " inner CRT slots from cylotomic poly " << hint_crt_values[0][0].size() << std::endl;

  // CRT interpolates to get Hint * LWE secrets as balanced mod-t values, where
  // t is the product of plaintext moduli. Then convert them wrt LWE modulus.
  BigInteger p = 1;
  for (auto pi : plaintext_moduli) {
    p *= rlwe::ConvertToBigInteger<RlweInteger, BigInteger>(pi->Modulus());
  }
  BigInteger p_half = p / 2;
  BigInteger lwe_modulus = BigInteger(SOL_NTT_PRIME);
  if (SOL_NTT_PRIME == 0){
    lwe_modulus = BigInteger(1) << params_.lwe_modulus_bit_size;
  }
  //std::cout << "THIS THIS LWE modulus " << lwe_modulus << std::endl;
  //BigInteger lwe_modulus = BigInteger(1) << params_.lwe_modulus_bit_size;
  std::vector<BigInteger> p_hats =
      rlwe::RnsModulusComplements<RlweModularInt, BigInteger>(plaintext_moduli);
  RLWE_ASSIGN_OR_RETURN(
      std::vector<RlweModularInt> p_hat_invs,
      crt_context_.MainPrimeModulusCrtFactors(num_linpir_plaintext_moduli - 1));

  std::vector<lwe::Vector> hint_vectors;
  hint_vectors.reserve(num_shards);
  for (int j = 0; j < num_shards; ++j) {
    RLWE_ASSIGN_OR_RETURN(
        std::vector<BigInteger> hint_values,
        (rlwe::CrtInterpolation<RlweModularInt, BigInteger>(
            hint_crt_values[j], plaintext_moduli, p_hats, p_hat_invs)));

    lwe::Vector hint = lwe::Vector::Zero(hint_values.size());
    for (int i = 0; i < hint_values.size(); ++i) {
      BigInteger x = hint_values[i] % p;
      hint[i] = static_cast<RlweInteger>(ConvertModulus(x, p, lwe_modulus, p_half));
    }
    hint_vectors.push_back(std::move(hint));
  }

  return hint_vectors;
}

}  // namespace hintless_simplepir
}  // namespace hintless_pir
