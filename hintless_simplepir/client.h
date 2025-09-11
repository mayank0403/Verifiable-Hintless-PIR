/*
 * Modified by Mayank Rathee
 * Copyright 2024 Google LLC.
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef HINTLESS_PIR_HINTLESS_SIMPLEPIR_CLIENT_H_
#define HINTLESS_PIR_HINTLESS_SIMPLEPIR_CLIENT_H_

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <tuple>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "hintless_simplepir/parameters.h"
#include "hintless_simplepir/serialization.pb.h"
#include "linpir/client.h"
#include "lwe/types.h"
#include "lwe/matrix_conversion.h"
#include "lwe/lwe_symmetric_encryption.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"

namespace hintless_pir {
namespace hintless_simplepir {

// The client part of the HintlessPir protocol.
class Client {
 public:
  // Creates a client from the given protocol parameters `params` and the
  // server's public parameters `public_params`.
  static absl::StatusOr<std::unique_ptr<Client>> Create(
      const Parameters& params,
      const HintlessPirServerPublicParams& public_params);

  // Returns the request for accessing database[index]. 
  absl::StatusOr<HintlessPirRequest> GenerateRequest(int64_t index);

  // Additionally returns secret key for VSPIR
  absl::StatusOr<std::pair<HintlessPirRequest, lwe::Vector>> GenerateRequest(int64_t index, lwe::Vector& lwe_sk_vsp_out);

  absl::StatusOr<std::pair<HintlessPirRequest, lwe::Vector>> GenerateRequestGivenAs(int64_t index, lwe::Vector& As, lwe::SymmetricLweKey& s_lwe);

  absl::StatusOr<std::pair<HintlessPirRequest, lwe::Vector>> GenerateRequestBSGS(int64_t index, lwe::Vector& lwe_sk_vsp_out);

  absl::StatusOr<std::pair<HintlessPirRequest, lwe::Vector>> GenerateRequestBSGSGivenAs(int64_t index, lwe::Vector& As, lwe::SymmetricLweKey& s_lwe);

  absl::StatusOr<HintlessPirRequest> PrepareLinPirGivenS_BSGS(lwe::SymmetricLweKey& s_lwe);

  absl::StatusOr<HintlessPirRequest> PrepareLinPirGivenS(lwe::SymmetricLweKey& s_lwe);

  absl::StatusOr<std::pair<HintlessPirRequest, lwe::Vector>> GenerateRequestGivenAsSkipLinPir(int64_t index, lwe::Vector& As, lwe::SymmetricLweKey& s_lwe);

  absl::StatusOr<std::tuple<VMatrix, lwe::Vector, lwe::SymmetricLweKey>> Compute_A_times_s();

  // Also returns sk to be passed to VSPIR and doesn't generate a SimplePIR query b/c that is handled by VSPIR
  absl::StatusOr<HintlessPirRequest> GenerateRequestWithoutQuery(int64_t index, lwe::Vector& lwe_sk_vsp_out, lwe::Vector& query_vsp_out);

  // Returns the retrieved record from the server response.
  absl::StatusOr<std::string> RecoverRecord(
      const HintlessPirResponse& response);

  absl::StatusOr<std::tuple<std::string, std::vector<lwe::Vector>, std::vector<std::vector<lwe::Integer>>>>  RecoverRecordAndGiveHsAndV(
    const HintlessPirResponse& response);

  absl::StatusOr<std::vector<lwe::Vector>> RecoverHsPreparePhase(
    const HintlessPirResponse& response);

  absl::StatusOr<std::tuple<std::string, std::vector<std::vector<lwe::Integer>>>>  RecoverRecordAndGiveV(
    const HintlessPirResponse& response, const std::vector<lwe::Vector>& decryption_parts);

  absl::StatusOr<std::tuple<std::string, lwe::Vector, std::vector<lwe::Integer>>> RecoverRecordAndGiveHsAndVNaive(
    const HintlessPirResponse& response, std::vector<lwe::Vector>& decryption_parts);

  absl::StatusOr<std::pair<std::string, lwe::Vector>> RecoverRecordGivenSimplePIRResp(
      const HintlessPirResponse& response, std::vector<std::vector<lwe::Integer>> res_in);

 private:
  using RlweInteger = Parameters::RlweInteger;
  using RlweModularInt = rlwe::MontgomeryInt<RlweInteger>;
  using RlweRnsContext = rlwe::RnsContext<RlweModularInt>;
  using RlwePrimeModulus = rlwe::PrimeModulus<RlweModularInt>;
  using LinPirClient = linpir::Client<RlweInteger>;

  // HintlessPir client state, which is cached for each request until the
  // corresponding response is received:
  //
  // 1) as in SimplePIR, a pair of indices (row_idx, col_idx) representing the
  // client's desired query index i = (row_idx * cols) + col_idx.
  //
  // 2) a PRNG seed expanding to the LinPir secret key for encrypting the LWE
  // secret used by the request.
  struct ClientState {
    int64_t row_idx;
    int64_t col_idx;
    std::string prng_seed_linpir_sk;
  };
  struct ClientStateNew {
    int64_t row_idx_begin;
    int64_t row_idx_end;
    int64_t col_idx;
    std::string prng_seed_linpir_sk;
  };

  explicit Client(
      Parameters params, absl::string_view prng_seed_lwe_query_pad,
      std::vector<std::unique_ptr<const RlweRnsContext>> rlwe_contexts,
      std::vector<const RlwePrimeModulus*> rlwe_moduli,
      std::vector<std::unique_ptr<LinPirClient>> linpir_clients,
      RlweRnsContext crt_context)
      : params_(std::move(params)),
        prng_seed_lwe_query_pad_(std::string(prng_seed_lwe_query_pad)),
        rlwe_contexts_(std::move(rlwe_contexts)),
        rlwe_moduli_(std::move(rlwe_moduli)),
        linpir_clients_(std::move(linpir_clients)),
        crt_context_(std::move(crt_context)) {}

  static std::vector<RlweInteger> EncodeLweVector(const lwe::Vector& lwe_vector,
                                                  RlweInteger lwe_modulus,
                                                  RlweInteger encode_modulus);

  // Encrypts the LWE secret vector using LinPir clients and update `request`
  // with the LinPir requests.
  absl::Status GenerateLinPirRequestInPlace(
      HintlessPirRequest& request, const lwe::Vector& lwe_secret) const;
    
  absl::Status GenerateLinPirRequestInPlaceBSGS(
    HintlessPirRequest& request, const lwe::Vector& lwe_secret) const;

  // CRT interpolates the LinPir responses to recover the LWE decryption parts,
  // which are the inner products hint * LWE secrets.
  absl::StatusOr<std::vector<lwe::Vector>> RecoverLweDecryptionParts(
      const HintlessPirResponse& response) const;

  const Parameters params_;

  // PRNG seed for generating the "A" matrix for LWE query ciphertext.
  std::string prng_seed_lwe_query_pad_;

  const std::vector<std::unique_ptr<const RlweRnsContext>> rlwe_contexts_;

  // The RLWE RNS moduli common to all LinPir clients.
  const std::vector<const RlwePrimeModulus*> rlwe_moduli_;

  const std::vector<std::unique_ptr<LinPirClient>> linpir_clients_;

  const RlweRnsContext crt_context_;

  // Per request state.
  ClientState state_;
  ClientStateNew state_new_;
};

}  // namespace hintless_simplepir
}  // namespace hintless_pir

#endif  // HINTLESS_PIR_HINTLESS_SIMPLEPIR_CLIENT_H_
