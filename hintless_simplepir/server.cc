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

#include "hintless_simplepir/server.h"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Eigen/Core"
#include "absl/memory/memory.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"

#include "hintless_simplepir/database_hwy.h"

#include "hintless_simplepir/parameters.h"
#include "hintless_simplepir/serialization.pb.h"
#include "hintless_simplepir/utils.h"
#include "lwe/lwe_symmetric_encryption.h"
#include "lwe/types.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/status_macros.h"

namespace hintless_pir {
namespace hintless_simplepir {

namespace {

// Returns an error if `params` uses an invalid PRNG type.
inline absl::Status CheckForValidPrngType(const Parameters& params) {
  if (!(params.prng_type == rlwe::PRNG_TYPE_HKDF ||
        params.prng_type == rlwe::PRNG_TYPE_CHACHA)) {
    return absl::InvalidArgumentError("Invalid PRNG type in `params`.");
  }
  return absl::OkStatus();
}

}  // namespace

absl::StatusOr<std::unique_ptr<Server>> Server::Create(
    const Parameters& params) {
  RLWE_RETURN_IF_ERROR(CheckForValidPrngType(params));

  // Create RLWE contexts, one per plaintext modulus in `ts`.
  auto const& rlwe_params = params.linpir_params;
  int num_linpir_instances = rlwe_params.ts.size();
  std::vector<std::unique_ptr<const RlweRnsContext>> rlwe_contexts;
  rlwe_contexts.reserve(num_linpir_instances);
  for (int i = 0; i < num_linpir_instances; ++i) {
    RLWE_ASSIGN_OR_RETURN(
        auto rlwe_context,
        RlweRnsContext::CreateForBfvFiniteFieldEncoding(
            rlwe_params.log_n, rlwe_params.qs, /*ps=*/{}, rlwe_params.ts[i]));
    rlwe_contexts.push_back(
        std::make_unique<const RlweRnsContext>(std::move(rlwe_context)));
  }

  // Create a Database object holding the database and hint matrices.
  RLWE_ASSIGN_OR_RETURN(auto database, Database::Create(params));

  return absl::WrapUnique(
      new Server(params, std::move(database), std::move(rlwe_contexts)));
}

absl::StatusOr<std::unique_ptr<Server>> Server::CreateWithRandomDatabaseRecordsFake(
    const Parameters& params) {
  RLWE_RETURN_IF_ERROR(CheckForValidPrngType(params));

  // Create RLWE contexts, one per plaintext modulus in `ts`.
  auto const& rlwe_params = params.linpir_params;
  int num_linpir_instances = rlwe_params.ts.size();
  std::vector<std::unique_ptr<const RlweRnsContext>> rlwe_contexts;
  rlwe_contexts.reserve(num_linpir_instances);
  for (int i = 0; i < num_linpir_instances; ++i) {
    //std::cout << i << " log n " << rlwe_params.log_n << " q " << rlwe_params.qs[0] << ", " << rlwe_params.qs[1] << " ts " << rlwe_params.ts[i] << std::endl;
    RLWE_ASSIGN_OR_RETURN(
        auto rlwe_context,
        RlweRnsContext::CreateForBfvFiniteFieldEncoding(
            rlwe_params.log_n, rlwe_params.qs, /*ps=*/{}, rlwe_params.ts[i]));
    rlwe_contexts.push_back(
        std::make_unique<const RlweRnsContext>(std::move(rlwe_context)));
  }
  std::cout << "Created RLWE contexts" << std::endl;
  // Create a Database holding random records.
  RLWE_ASSIGN_OR_RETURN(auto database, Database::CreateRandomFake(params));

  return absl::WrapUnique(
      new Server(params, std::move(database), std::move(rlwe_contexts)));
}

absl::StatusOr<std::unique_ptr<Server>> Server::CreateWithRandomDatabaseRecords(
    const Parameters& params) {
  RLWE_RETURN_IF_ERROR(CheckForValidPrngType(params));

  // Create RLWE contexts, one per plaintext modulus in `ts`.
  auto const& rlwe_params = params.linpir_params;
  int num_linpir_instances = rlwe_params.ts.size();
  std::vector<std::unique_ptr<const RlweRnsContext>> rlwe_contexts;
  rlwe_contexts.reserve(num_linpir_instances);
  for (int i = 0; i < num_linpir_instances; ++i) {
    //std::cout << i << " log n " << rlwe_params.log_n << " q " << rlwe_params.qs[0] << ", " << rlwe_params.qs[1] << " ts " << rlwe_params.ts[i] << std::endl;
    RLWE_ASSIGN_OR_RETURN(
        auto rlwe_context,
        RlweRnsContext::CreateForBfvFiniteFieldEncoding(
            rlwe_params.log_n, rlwe_params.qs, /*ps=*/{}, rlwe_params.ts[i]));
    rlwe_contexts.push_back(
        std::make_unique<const RlweRnsContext>(std::move(rlwe_context)));
  }
  std::cout << "Created RLWE contexts" << std::endl;
  // Create a Database holding random records.
  RLWE_ASSIGN_OR_RETURN(auto database, Database::CreateRandom(params));

  return absl::WrapUnique(
      new Server(params, std::move(database), std::move(rlwe_contexts)));
}

absl::StatusOr<std::unique_ptr<Server>> Server::CreateWithFastDatabaseRecords(
    const Parameters& params) {
  RLWE_RETURN_IF_ERROR(CheckForValidPrngType(params));

  // Create RLWE contexts, one per plaintext modulus in `ts`.
  auto const& rlwe_params = params.linpir_params;
  int num_linpir_instances = rlwe_params.ts.size();
  std::vector<std::unique_ptr<const RlweRnsContext>> rlwe_contexts;
  rlwe_contexts.reserve(num_linpir_instances);
  for (int i = 0; i < num_linpir_instances; ++i) {
    //std::cout << i << std::endl;
    RLWE_ASSIGN_OR_RETURN(
        auto rlwe_context,
        RlweRnsContext::CreateForBfvFiniteFieldEncoding(
            rlwe_params.log_n, rlwe_params.qs, /*ps=*/{}, rlwe_params.ts[i]));
    rlwe_contexts.push_back(
        std::make_unique<const RlweRnsContext>(std::move(rlwe_context)));
  }

  // Create a Databas holding random records.
  RLWE_ASSIGN_OR_RETURN(auto database, Database::CreateFast(params));

  return absl::WrapUnique(
      new Server(params, std::move(database), std::move(rlwe_contexts)));
}

absl::Status Server::GeneratePublicParams() {
  int num_linpir_instances = params_.linpir_params.ts.size();
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    // Sample PRNG seeds for LWE "A" matrix and LinPIR.
    RLWE_ASSIGN_OR_RETURN(prng_seed_lwe_query_pad_,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    prng_seed_linpir_ct_pads_.clear();
    prng_seed_linpir_ct_pads_.resize(num_linpir_instances);
    for (int i = 0; i < num_linpir_instances; ++i) {
      RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_ct_pads_[i],
                            rlwe::SingleThreadHkdfPrng::GenerateSeed());
    }
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_gk_pad_,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_gk_pad_baby_,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_gk_pad_giant_,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    
    // Generate the LWE "A" matrix.
    RLWE_ASSIGN_OR_RETURN(auto prng, rlwe::SingleThreadHkdfPrng::Create(
                                         prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        auto pad,
        lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim, prng.get()));
    lwe_query_pad_ = std::make_unique<const lwe::Matrix>(std::move(pad));
  } else {
    RLWE_ASSIGN_OR_RETURN(prng_seed_lwe_query_pad_,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    prng_seed_linpir_ct_pads_.clear();
    prng_seed_linpir_ct_pads_.resize(num_linpir_instances);
    for (int i = 0; i < num_linpir_instances; ++i) {
      RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_ct_pads_[i],
                            rlwe::SingleThreadChaChaPrng::GenerateSeed());
    }
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_gk_pad_,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_gk_pad_baby_,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_seed_linpir_gk_pad_giant_,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());

    RLWE_ASSIGN_OR_RETURN(auto prng, rlwe::SingleThreadChaChaPrng::Create(
                                         prng_seed_lwe_query_pad_));
    RLWE_ASSIGN_OR_RETURN(
        auto pad,
        lwe::ExpandPad(params_.db_cols, params_.lwe_secret_dim, prng.get()));

    /*
    // if solinas_32 is true, take modulo of pad
    // TODO: Note that this doesn't guarantee uniform pad over Z_q; need rejection sampling for that; for this prototype, we are ignoring this.
    // Sampling of A doesn't affect runtime costs, so ignoring rejection sampling of A for now.
    uint64_t rows = pad.rows();
    uint64_t cols = pad.cols();
    for (uint64_t i = 0; i < rows; i++) {
      for (uint64_t j = 0; j < cols; j++) {
        pad(i, j) = fast_mod_max_64(pad(i, j));
      }
    }
    */

    lwe_query_pad_ = std::make_unique<const lwe::Matrix>(std::move(pad));
  }
  return absl::OkStatus();
}

namespace {

// Given `matrix` with mod-q entries, returns `matrix` mod p, where modular
// numbers are in balanced representation.

template <typename Integer>
std::vector<std::vector<Integer>> EncodeLweMatrix(
    const Database::LweMatrix& matrix, Integer q, Integer p) {
  Integer q_half = q >> 1;
  int num_rows = matrix.size();
  int num_cols = matrix[0].size();
  std::vector<std::vector<Integer>> matrix_mod_p(num_rows);
  for (int i = 0; i < num_rows; ++i) {
    matrix_mod_p[i].reserve(num_cols);
    for (int j = 0; j < num_cols; ++j) {
      Integer x = static_cast<Integer>(matrix[i][j]);
      matrix_mod_p[i].push_back(ConvertModulus(x, q, p, q_half));
    }
  }
  return matrix_mod_p;
}


}  // namespace

absl::Status Server::Preprocess() {
  // Refresh the PRNG seeds.
  RLWE_RETURN_IF_ERROR(GeneratePublicParams()); // Generates gk_pad_seed for BSGS and also non-BSGS case

  std::cout << "Generated public params" << std::endl;

  // Make sure the hint is up to date.
  RLWE_RETURN_IF_ERROR(database_->UpdateLweQueryPad(lwe_query_pad_.get()));
#ifdef FAKE_RUN
  RLWE_RETURN_IF_ERROR(database_->UpdateHintsFake());
#else
  RLWE_RETURN_IF_ERROR(database_->UpdateHints());
#endif
  //std::cout << "Updated hints" << std::endl;

  RlweInteger lwe_modulus = RlweInteger{1} << params_.lwe_modulus_bit_size;
  if (SOL_NTT_PRIME > 0) {
    lwe_modulus = SOL_NTT_PRIME;
  }
  size_t num_shards = database_->NumShards();

  // Create LinPir databases (holding the preprocessed hints) and servers.
  linpir_servers_.clear();
  linpir_servers_.resize(rlwe_contexts_.size());
  linpir_databases_.clear();
  linpir_databases_.resize(rlwe_contexts_.size());
  for (int k = 0; k < rlwe_contexts_.size(); ++k) {
    std::cout << "Creating LinPIR databases and servers" << std::endl;
    RlweInteger plaintext_modulus = rlwe_contexts_[k]->PlaintextModulus();

    // One LinPir database per shard, for the current plaintext modulus.
    std::vector<std::unique_ptr<LinPirDatabase>> linpir_databases_mod_tk;
    linpir_databases_mod_tk.reserve(num_shards);


    for (const Database::LweMatrix& hint : database_->Hints()) {
    
      std::vector<std::vector<RlweInteger>> hint_mod_tk =
          EncodeLweMatrix(hint, lwe_modulus, plaintext_modulus);
      RLWE_ASSIGN_OR_RETURN(
          auto linpir_database,
          LinPirDatabase::Create(params_.linpir_params, rlwe_contexts_[k].get(),
                                 hint_mod_tk));
      linpir_databases_mod_tk.push_back(std::move(linpir_database));
    }
    std::cout << "Created LinPIR database" << std::endl;
    std::vector<LinPirDatabase*> linpir_databases_ptrs;
    std::transform(linpir_databases_mod_tk.begin(),
                   linpir_databases_mod_tk.end(),
                   std::back_inserter(linpir_databases_ptrs),
                   [](auto& ptr) { return ptr.get(); });
    //std::cout << "Transform done" << std::endl;
    RLWE_ASSIGN_OR_RETURN(
        auto linpir_server_mod_tk,
        LinPirServer::Create(params_.linpir_params, rlwe_contexts_[k].get(),
                             linpir_databases_ptrs,
                             prng_seed_linpir_ct_pads_[k],
                             prng_seed_linpir_gk_pad_));

    RLWE_RETURN_IF_ERROR(linpir_server_mod_tk->Preprocess());

    linpir_databases_[k] = std::move(linpir_databases_mod_tk);
    linpir_servers_[k] = std::move(linpir_server_mod_tk);
  }

  return absl::OkStatus();
}

absl::Status Server::PreprocessBSGS() {
  // Refresh the PRNG seeds.
  RLWE_RETURN_IF_ERROR(GeneratePublicParams()); // Generates gk_pad_seed for BSGS and also non-BSGS case
  
  //std::cout << "Generated public params" << std::endl;

  // Make sure the hint is up to date.
  RLWE_RETURN_IF_ERROR(database_->UpdateLweQueryPad(lwe_query_pad_.get()));
#ifdef FAKE_RUN
  RLWE_RETURN_IF_ERROR(database_->UpdateHintsFake());
#else
  RLWE_RETURN_IF_ERROR(database_->UpdateHints());
#endif
  //RLWE_RETURN_IF_ERROR(database_->UpdateHints());
  //std::cout << "Updated hints" << std::endl;

  RlweInteger lwe_modulus = RlweInteger{1} << params_.lwe_modulus_bit_size;
  if (SOL_NTT_PRIME > 0) {
    lwe_modulus = SOL_NTT_PRIME;
  }
  size_t num_shards = database_->NumShards();

  // Create LinPir databases (holding the preprocessed hints) and servers.
  linpir_servers_.clear();
  linpir_servers_.resize(rlwe_contexts_.size());
  linpir_databases_.clear();
  linpir_databases_.resize(rlwe_contexts_.size());
  for (int k = 0; k < rlwe_contexts_.size(); ++k) {
    //std::cout << "Creating LinPIR databases and servers" << std::endl;
    RlweInteger plaintext_modulus = rlwe_contexts_[k]->PlaintextModulus();

    // One LinPir database per shard, for the current plaintext modulus.
    std::vector<std::unique_ptr<LinPirDatabase>> linpir_databases_mod_tk;
    linpir_databases_mod_tk.reserve(num_shards);


    for (const Database::LweMatrix& hint : database_->Hints()) {
    
      std::vector<std::vector<RlweInteger>> hint_mod_tk =
          EncodeLweMatrix(hint, lwe_modulus, plaintext_modulus);
      RLWE_ASSIGN_OR_RETURN(
          auto linpir_database,
          LinPirDatabase::CreateBSGS(params_.linpir_params, rlwe_contexts_[k].get(),
                                 hint_mod_tk));
      linpir_databases_mod_tk.push_back(std::move(linpir_database));
    }
    //std::cout << "Created LinPIR database" << std::endl;
    std::vector<LinPirDatabase*> linpir_databases_ptrs;
    std::transform(linpir_databases_mod_tk.begin(),
                   linpir_databases_mod_tk.end(),
                   std::back_inserter(linpir_databases_ptrs),
                   [](auto& ptr) { return ptr.get(); });
    //std::cout << "Transform done" << std::endl;
    RLWE_ASSIGN_OR_RETURN(
        auto linpir_server_mod_tk,
        LinPirServer::CreateBSGS(params_.linpir_params, rlwe_contexts_[k].get(),
                             linpir_databases_ptrs,
                             prng_seed_linpir_ct_pads_[k],
                             prng_seed_linpir_gk_pad_baby_,
                             prng_seed_linpir_gk_pad_giant_));
    
    RLWE_RETURN_IF_ERROR(linpir_server_mod_tk->PreprocessBSGS());

    linpir_databases_[k] = std::move(linpir_databases_mod_tk);
    linpir_servers_[k] = std::move(linpir_server_mod_tk);
  }

  return absl::OkStatus();
}

absl::StatusOr<HintlessPirResponse> Server::HandleRequest(
    const HintlessPirRequest& request, bool skip_pir_query) {
  if (!IsPreprocessed()) {
    return absl::FailedPreconditionError("Server has not been preprocessed.");
  }

  HintlessPirResponse response;
  double start_lwe, start_Rlwe, end_lwe, end_Rlwe;
  start_lwe = currentDateTime();
  // Handle the LWE part of the request.
  if (!skip_pir_query){
    std::cout<<"Deser LWE ciphertext"<<std::endl;
    Database::LweVector ct_query_vector =
        DeserializeLweCiphertext(request.ct_query_vector());
    RLWE_ASSIGN_OR_RETURN(std::vector<Database::LweVector> ct_records,
                          database_->InnerProductWith(ct_query_vector));
    
    //std::cout<<"Ser LWE ciphertext"<<std::endl;
    for (auto& ct_record : ct_records) {
      *response.add_ct_records() = SerializeLweCiphertext(ct_record);
    }
  }
  end_lwe = currentDateTime();
  std::cout << "[==> TIMER  <==] Server-only Online time for D*u (LWE): " << (end_lwe-start_lwe) << " ms | " << (end_lwe-start_lwe)/1000 << " sec" << std::endl;
  
  start_Rlwe = currentDateTime();
  // Handle the LinPIR requests.
  int num_linpir_requests = request.linpir_ct_bs_size();
  if (num_linpir_requests != linpir_servers_.size()) {
    return absl::InvalidArgumentError(
        "`request` contains unexpected number of LinPir requests.");
  }
  for (int k = 0; k < num_linpir_requests; ++k) {
    std::cout<<"LinPir handle request"<<std::endl;
    RLWE_ASSIGN_OR_RETURN(LinPirResponse linpir_response,
                          linpir_servers_[k]->HandleRequest(
                              request.linpir_ct_bs(k), request.linpir_gk_bs()));
    *response.add_linpir_responses() = std::move(linpir_response);
  }
  end_Rlwe = currentDateTime();
  std::cout << "[==> TIMER  <==] Server-only Online time for H*s (R-LWE): " << (end_Rlwe-start_Rlwe) << " ms | " << (end_Rlwe-start_Rlwe)/1000 << " sec" << std::endl;

  return response;
}

absl::StatusOr<HintlessPirResponse> Server::HandleRequestBSGS(
    const HintlessPirRequest& request, bool skip_pir_query) {
  if (!IsPreprocessed()) {
    return absl::FailedPreconditionError("Server has not been preprocessed.");
  }

  HintlessPirResponse response;
  double start_lwe, start_Rlwe, end_lwe, end_Rlwe;
  start_lwe = currentDateTime();
  // Handle the LWE part of the request.
  if (!skip_pir_query){
    //std::cout<<"Deser LWE ciphertext"<<std::endl;
    Database::LweVector ct_query_vector =
        DeserializeLweCiphertext(request.ct_query_vector());
    RLWE_ASSIGN_OR_RETURN(std::vector<Database::LweVector> ct_records,
                          database_->InnerProductWith(ct_query_vector));
    
    //std::cout<<"Ser LWE ciphertext"<<std::endl;
    for (auto& ct_record : ct_records) {
      *response.add_ct_records() = SerializeLweCiphertext(ct_record);
    }
  }
  end_lwe = currentDateTime();
  std::cout << BLUE << "[==> TIMER  <==] Server-only Online time for D*u (LWE): " << (end_lwe-start_lwe) << " ms | " << (end_lwe-start_lwe)/1000 << " sec" << END << std::endl;
  
  start_Rlwe = currentDateTime();
  // Handle the LinPIR requests.
  int num_linpir_requests = request.linpir_ct_bs_size();
  if (num_linpir_requests != linpir_servers_.size()) {
    return absl::InvalidArgumentError(
        "`request` contains unexpected number of LinPir requests.");
  }
  for (int k = 0; k < num_linpir_requests; ++k) {
    //std::cout<<"LinPir handle request"<<std::endl;
    RLWE_ASSIGN_OR_RETURN(LinPirResponse linpir_response,
                          linpir_servers_[k]->HandleRequestBSGS(
                              request.linpir_ct_bs(k), request.linpir_gk_bs_baby(), request.linpir_gk_bs_giant()));
    *response.add_linpir_responses() = std::move(linpir_response);
  }
  end_Rlwe = currentDateTime();
  std::cout << BLUE << "[==> TIMER  <==] Server-only Online time for H*s (R-LWE): " << (end_Rlwe-start_Rlwe) << " ms | " << (end_Rlwe-start_Rlwe)/1000 << " sec" << END << std::endl;

  return response;
}

absl::StatusOr<HintlessPirResponse> Server::HandlePrepareRequest(
    const HintlessPirRequest& request) {
  if (!IsPreprocessed()) {
    return absl::FailedPreconditionError("Server has not been preprocessed.");
  }
  HintlessPirResponse response;

  // Handle the LinPIR requests.
  int num_linpir_requests = request.linpir_ct_bs_size();
  if (num_linpir_requests != linpir_servers_.size()) {
    return absl::InvalidArgumentError(
        "`request` contains unexpected number of LinPir requests.");
  }
  for (int k = 0; k < num_linpir_requests; ++k) {
    //std::cout<<"LinPir handle request"<<std::endl;
    RLWE_ASSIGN_OR_RETURN(LinPirResponse linpir_response,
                          linpir_servers_[k]->HandleRequest(
                              request.linpir_ct_bs(k), request.linpir_gk_bs()));
    *response.add_linpir_responses() = std::move(linpir_response);
  }

  return response;
}

absl::StatusOr<HintlessPirResponse> Server::HandlePrepareRequestBSGS(
    const HintlessPirRequest& request) {
  if (!IsPreprocessed()) {
    return absl::FailedPreconditionError("Server has not been preprocessed.");
  }
  HintlessPirResponse response;

  // Handle the LinPIR requests.
  int num_linpir_requests = request.linpir_ct_bs_size();
  if (num_linpir_requests != linpir_servers_.size()) {
    return absl::InvalidArgumentError(
        "`request` contains unexpected number of LinPir requests.");
  }
  for (int k = 0; k < num_linpir_requests; ++k) {
    //std::cout<<"LinPir handle request"<<std::endl;
    RLWE_ASSIGN_OR_RETURN(LinPirResponse linpir_response,
                          linpir_servers_[k]->HandleRequestBSGS(
                              request.linpir_ct_bs(k), request.linpir_gk_bs_baby(), request.linpir_gk_bs_giant()));
    *response.add_linpir_responses() = std::move(linpir_response);
  }

  return response;
}

absl::StatusOr<HintlessPirResponse> Server::HandleRequestSkipLinPir(
    const HintlessPirRequest& request) {
  if (!IsPreprocessed()) {
    return absl::FailedPreconditionError("Server has not been preprocessed.");
  }

  HintlessPirResponse response;
  // Handle the LWE part of the request.
  
  //std::cout<<"Deser LWE ciphertext"<<std::endl;
  Database::LweVector ct_query_vector =
      DeserializeLweCiphertext(request.ct_query_vector());
  RLWE_ASSIGN_OR_RETURN(std::vector<Database::LweVector> ct_records,
                        database_->InnerProductWith(ct_query_vector));
  
  //std::cout<<"Ser LWE ciphertext"<<std::endl;
  for (auto& ct_record : ct_records) {
    *response.add_ct_records() = SerializeLweCiphertext(ct_record);
  }

  return response;
}

HintlessPirServerPublicParams Server::GetPublicParams() const {
  HintlessPirServerPublicParams output;
  output.set_prng_seed_lwe_query_pad(prng_seed_lwe_query_pad_);
  for (auto const& prng_seed : prng_seed_linpir_ct_pads_) {
    *output.add_prng_seed_linpir_ct_pads() = prng_seed;
  }
  output.set_prng_seed_linpir_gk_pad(prng_seed_linpir_gk_pad_);
  output.set_prng_seed_linpir_gk_pad_baby(prng_seed_linpir_gk_pad_baby_);
  output.set_prng_seed_linpir_gk_pad_giant(prng_seed_linpir_gk_pad_giant_);

  return output;
}

absl::StatusOr<std::vector<Database::LweMatrix>> Server::TransposeProduct(
    const Database::LweMatrix& A, bool secondary) {
  
  return database_->TransposeMatrixProduct(A, secondary);
}

absl::StatusOr<std::vector<Database::LweMatrix64>> Server::TransposeProduct64(
    const Database::LweMatrix64& A, bool secondary) {
  
  return database_->TransposeMatrixProduct64(A, secondary);
}

absl::StatusOr<std::vector<Database::LweMatrix>> Server::DirectProduct(
    const Database::LweMatrix& A, bool secondary) {
  
  return database_->DirectMatrixProduct(A, secondary);
}

std::vector<Database::LweMatrix> Server::TransposeProductFake() {
  
  return database_->TransposeMatrixProductFake();
}

std::vector<Database::LweMatrix64> Server::TransposeProductFake64() {
  
  return database_->TransposeMatrixProductFake64();
}

std::vector<Database::LweMatrix> Server::DirectProductFake() {
  
  return database_->DirectMatrixProductFake();
}

// Outermost layer is shards, then D^T * each ciphertext
absl::StatusOr<std::vector<std::vector<Database::LweVector>>> Server::TransposeProductCiphertexts(
    const std::vector<Database::LweVector>& cts) {
  
  return database_->TransposeMatrixCiphertextProducts(cts);
}

// Outermost layer is shards, then D^T * each ciphertext
absl::StatusOr<std::vector<std::vector<Database::LweVector64>>> Server::TransposeProductCiphertexts64(
    const std::vector<Database::LweVector64>& cts) {
  
  return database_->TransposeMatrixCiphertextProducts64(cts);
}

}  // namespace hintless_simplepir
}  // namespace hintless_pir
