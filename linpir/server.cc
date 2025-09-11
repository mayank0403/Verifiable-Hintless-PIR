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

#include "linpir/server.h"

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/memory/memory.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "google/protobuf/repeated_ptr_field.h"
#include "linpir/database.h"
#include "linpir/parameters.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/status_macros.h"

int small_modexp_2n(int base, int exp, int mod) {
  int result = 1;
  for (int i = 0; i < exp; i++) {
    result = (result * base) % (2 * mod); // 2 * mod because in X^n + 1 cylotomic polynomial, we have that X^{2*n} = 1 unlike X^n = -1.
  }
  return result;
}

namespace hintless_pir {
namespace linpir {

template <typename RlweInteger>
absl::StatusOr<std::unique_ptr<Server<RlweInteger>>>
Server<RlweInteger>::Create(
    const RlweParameters<RlweInteger>& parameters,
    const RnsContext* rns_context,
    const std::vector<Database<RlweInteger>*>& databases,
    absl::string_view prng_seed_ct_pad, absl::string_view prng_seed_gk_pad) {
  if (!(parameters.prng_type == rlwe::PRNG_TYPE_HKDF ||
        parameters.prng_type == rlwe::PRNG_TYPE_CHACHA)) {
    return absl::InvalidArgumentError("Invalid `prng_type`.");
  }
  if (rns_context == nullptr) {
    return absl::InvalidArgumentError("`rns_context` must not be null.");
  }

  auto rns_moduli = rns_context->MainPrimeModuli();
  int level = rns_moduli.size() - 1;
  RLWE_ASSIGN_OR_RETURN(auto q_hats,
                        rns_context->MainPrimeModulusComplements(level));
  RLWE_ASSIGN_OR_RETURN(auto q_hat_invs,
                        rns_context->MainPrimeModulusCrtFactors(level));
  RLWE_ASSIGN_OR_RETURN(
      RnsGadget rns_gadget,
      RnsGadget::Create(parameters.log_n, parameters.gadget_log_bs, q_hats,
                        q_hat_invs, rns_moduli));
  RLWE_ASSIGN_OR_RETURN(
      auto rns_error_params,
      RnsErrorParams::Create(
          parameters.log_n, rns_moduli, /*aux_moduli=*/{},
          std::log2(static_cast<double>(rns_context->PlaintextModulus())),
          std::sqrt(parameters.error_variance)));

  return absl::WrapUnique(new Server<RlweInteger>(
      parameters, std::string(prng_seed_ct_pad), std::string(prng_seed_gk_pad),
      rns_context, std::move(rns_moduli), std::move(rns_gadget),
      std::move(rns_error_params), databases));
}

template <typename RlweInteger>
absl::StatusOr<std::unique_ptr<Server<RlweInteger>>>
Server<RlweInteger>::CreateBSGS(
    const RlweParameters<RlweInteger>& parameters,
    const RnsContext* rns_context,
    const std::vector<Database<RlweInteger>*>& databases,
    absl::string_view prng_seed_ct_pad, absl::string_view prng_seed_gk_pad_baby, absl::string_view prng_seed_gk_pad_giant) {
  if (!(parameters.prng_type == rlwe::PRNG_TYPE_HKDF ||
        parameters.prng_type == rlwe::PRNG_TYPE_CHACHA)) {
    return absl::InvalidArgumentError("Invalid `prng_type`.");
  }
  if (rns_context == nullptr) {
    return absl::InvalidArgumentError("`rns_context` must not be null.");
  }

  auto rns_moduli = rns_context->MainPrimeModuli();
  int level = rns_moduli.size() - 1;
  RLWE_ASSIGN_OR_RETURN(auto q_hats,
                        rns_context->MainPrimeModulusComplements(level));
  RLWE_ASSIGN_OR_RETURN(auto q_hat_invs,
                        rns_context->MainPrimeModulusCrtFactors(level));
  RLWE_ASSIGN_OR_RETURN(
      RnsGadget rns_gadget,
      RnsGadget::Create(parameters.log_n, parameters.gadget_log_bs, q_hats,
                        q_hat_invs, rns_moduli));
  RLWE_ASSIGN_OR_RETURN(
      auto rns_error_params,
      RnsErrorParams::Create(
          parameters.log_n, rns_moduli, /*aux_moduli=*/{},
          std::log2(static_cast<double>(rns_context->PlaintextModulus())),
          std::sqrt(parameters.error_variance)));

  return absl::WrapUnique(new Server<RlweInteger>(
      parameters, std::string(prng_seed_ct_pad), std::string(prng_seed_gk_pad_baby), std::string(prng_seed_gk_pad_giant),
      rns_context, std::move(rns_moduli), std::move(rns_gadget),
      std::move(rns_error_params), databases));
}

template <typename RlweInteger>
absl::StatusOr<std::unique_ptr<Server<RlweInteger>>>
Server<RlweInteger>::Create(
    const RlweParameters<RlweInteger>& parameters,
    const RnsContext* rns_context,
    const std::vector<Database<RlweInteger>*>& databases) {
  // Sample PRNG seeds for the query vector and the Galois key.
  std::string prng_seed_ct_pad, prng_seed_gk_pad;
  if (parameters.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_seed_ct_pad,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_seed_gk_pad,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
  } else if (parameters.prng_type == rlwe::PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(prng_seed_ct_pad,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_seed_gk_pad,
                          rlwe::SingleThreadChaChaPrng::GenerateSeed());
  } else {
    return absl::InvalidArgumentError("Invalid `prng_type`.");
  }

  return Server<RlweInteger>::Create(parameters, rns_context, databases,
                                     prng_seed_ct_pad, prng_seed_gk_pad);
}

template <typename RlweInteger>
absl::Status Server<RlweInteger>::Preprocess() {
  //std::cout << "Preprocessing LinPIR Server" << std::endl;
  ct_pads_.clear();
  ct_sub_pad_digits_.clear();
  gk_pads_.clear();

  // Create PRNGs.
  std::unique_ptr<rlwe::SecurePrng> prng_ct, prng_gk;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(
        prng_ct, rlwe::SingleThreadHkdfPrng::Create(prng_seed_ct_pad_));
    RLWE_ASSIGN_OR_RETURN(
        prng_gk, rlwe::SingleThreadHkdfPrng::Create(prng_seed_gk_pad_));
  } else {
    RLWE_ASSIGN_OR_RETURN(
        prng_ct, rlwe::SingleThreadChaChaPrng::Create(prng_seed_ct_pad_));
    RLWE_ASSIGN_OR_RETURN(
        prng_gk, rlwe::SingleThreadChaChaPrng::Create(prng_seed_gk_pad_));
  }

  // Expand seed to the "a" part of Enc(query vector)
  int log_n = rns_context_->LogN();
  int gadget_dim = rns_gadget_.Dimension();
  RLWE_ASSIGN_OR_RETURN(auto ct_pad, RnsPolynomial::SampleUniform(
                                         log_n, prng_ct.get(), rns_moduli_));
  RLWE_RETURN_IF_ERROR(ct_pad.NegateInPlace(rns_moduli_));

  // Create the "a" part of Galois key
  RLWE_ASSIGN_OR_RETURN(gk_pads_, RnsGaloisKey::SampleRandomPad(
                                      gadget_dim, log_n, rns_moduli_,
                                      prng_seed_gk_pad_, params_.prng_type));

  // Precompute the "a" part of Enc(s << i) and the digits used to generate
  // Enc(s << i).
  int num_rotations = params_.rows_per_block / 2;
  ct_pads_.reserve(num_rotations);
  ct_pads_.push_back(std::move(ct_pad));
  ct_sub_pad_digits_.reserve(num_rotations);

  int curr_power = 1;
  int cyclotomic_order = 1 << (log_n + 1);
  for (int i = 1; i < num_rotations; ++i) {
    curr_power = (curr_power * 5) % cyclotomic_order;
    // ct[i-1].a(X^5)
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial prev_sub_a,
                          ct_pads_[i - 1].Substitute(5, rns_moduli_));

    // g^-1(ct[i-1].a(X^5))
    if (prev_sub_a.IsNttForm()) {
      RLWE_RETURN_IF_ERROR(prev_sub_a.ConvertToCoeffForm(rns_moduli_));
    }
    RLWE_ASSIGN_OR_RETURN(auto prev_sub_a_digits,
                          rns_gadget_.Decompose(prev_sub_a, rns_moduli_));
    for (auto& digit : prev_sub_a_digits) {
      RLWE_RETURN_IF_ERROR(digit.ConvertToNttForm(rns_moduli_));
    }

    // g^-1(ct[i-1].a(X^5))^T * gk.a
    RLWE_ASSIGN_OR_RETURN(
        auto curr_a,
        RnsPolynomial::CreateZero(log_n, rns_moduli_, /*is_ntt=*/true));
    for (int j = 0; j < prev_sub_a_digits.size(); ++j) {
      RLWE_RETURN_IF_ERROR(curr_a.FusedMulAddInPlace(prev_sub_a_digits[j],
                                                     gk_pads_[j], rns_moduli_));
    }
    ct_pads_.push_back(std::move(curr_a));
    ct_sub_pad_digits_.push_back(std::move(prev_sub_a_digits));
  }

  // Preprocess the databases using the "a" part of Enc(s << i).
  for (auto const& database : databases_) {
    RLWE_RETURN_IF_ERROR(database->Preprocess(ct_pads_));
  }
  return absl::OkStatus();
}

template <typename RlweInteger>
absl::Status Server<RlweInteger>::PreprocessBSGS() {
  //std::cout << "Preprocessing LinPIR Server" << std::endl;
  ct_pads_baby_.clear();
  ct_sub_pad_baby_digits_.clear();
  ct_pads_giant_.clear();
  ct_sub_pad_giant_digits_.clear();
  gk_pads_baby_.clear();
  gk_pads_giant_.clear();

  // Create PRNGs.
  std::unique_ptr<rlwe::SecurePrng> prng_ct, prng_gk_baby, prng_gk_giant;
  if (params_.prng_type == rlwe::PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(
        prng_ct, rlwe::SingleThreadHkdfPrng::Create(prng_seed_ct_pad_));
    RLWE_ASSIGN_OR_RETURN(
        prng_gk_baby, rlwe::SingleThreadHkdfPrng::Create(prng_seed_gk_pad_baby_));
    RLWE_ASSIGN_OR_RETURN(
        prng_gk_giant, rlwe::SingleThreadHkdfPrng::Create(prng_seed_gk_pad_giant_));
  } else {
    RLWE_ASSIGN_OR_RETURN(
        prng_ct, rlwe::SingleThreadChaChaPrng::Create(prng_seed_ct_pad_));
    RLWE_ASSIGN_OR_RETURN(
        prng_gk_baby, rlwe::SingleThreadChaChaPrng::Create(prng_seed_gk_pad_baby_));
    RLWE_ASSIGN_OR_RETURN(
        prng_gk_giant, rlwe::SingleThreadChaChaPrng::Create(prng_seed_gk_pad_giant_));
  }

  // Expand seed to the "a" part of Enc(query vector)
  int log_n = rns_context_->LogN();
  int gadget_dim = rns_gadget_.Dimension();
  RLWE_ASSIGN_OR_RETURN(auto ct_pad, RnsPolynomial::SampleUniform(
                                         log_n, prng_ct.get(), rns_moduli_));
  RLWE_RETURN_IF_ERROR(ct_pad.NegateInPlace(rns_moduli_));

  // Create the "a" part of Galois keys
  RLWE_ASSIGN_OR_RETURN(gk_pads_baby_, RnsGaloisKey::SampleRandomPad(
                                      gadget_dim, log_n, rns_moduli_,
                                      prng_seed_gk_pad_baby_, params_.prng_type));
  RLWE_ASSIGN_OR_RETURN(gk_pads_giant_, RnsGaloisKey::SampleRandomPad(
                                      gadget_dim, log_n, rns_moduli_,
                                      prng_seed_gk_pad_giant_, params_.prng_type));

  // Precompute the "a" part of Enc(s << i) and the digits used to generate
  // Enc(s << i).
  
  int num_blocks = DivideAndRoundUp(params_.db_rows, params_.rows_per_block);
  int num_rotations = params_.rows_per_block / 2;
  auto baby_giant = find_optimal_bsgs(num_rotations, num_blocks);
  int baby_steps_total = baby_giant.first;
  int giant_steps_total = baby_giant.second;
  //int giant_steps_total = std::ceil(std::sqrt(num_rotations) / params_.baby_steps_bias_factor); // like row count
  //int baby_steps_total = std::ceil((float)num_rotations / giant_steps_total); // like column count
  int baby_step_power = 5;
  int giant_step_power = small_modexp_2n(baby_step_power, baby_steps_total, 1 << params_.log_n); 
  
  ct_pads_baby_.reserve(baby_steps_total);
  ct_pads_baby_.push_back(std::move(ct_pad));
  ct_sub_pad_baby_digits_.reserve(baby_steps_total);
  // TODO: DO SIMILAR FOR GIANT STEPS (I think size will be giant_steps_total * linpir_instances; confirm this)

  //int curr_power = 1;
  //int cyclotomic_order = 1 << (log_n + 1);
  // std::cout << "Performing baby step rotations" << std::endl;
  for (int i = 1; i < baby_steps_total; ++i) {
    //curr_power = (curr_power * 5) % cyclotomic_order;
    // ct[i-1].a(X^5)
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial prev_sub_a,
                          ct_pads_baby_[i - 1].Substitute(5, rns_moduli_));

    // g^-1(ct[i-1].a(X^5))
    if (prev_sub_a.IsNttForm()) {
      RLWE_RETURN_IF_ERROR(prev_sub_a.ConvertToCoeffForm(rns_moduli_));
    }
    RLWE_ASSIGN_OR_RETURN(auto prev_sub_a_digits,
                          rns_gadget_.Decompose(prev_sub_a, rns_moduli_));
    for (auto& digit : prev_sub_a_digits) {
      RLWE_RETURN_IF_ERROR(digit.ConvertToNttForm(rns_moduli_));
    }

    // g^-1(ct[i-1].a(X^5))^T * gk.a
    RLWE_ASSIGN_OR_RETURN(
        auto curr_a,
        RnsPolynomial::CreateZero(log_n, rns_moduli_, /*is_ntt=*/true));
    for (int j = 0; j < prev_sub_a_digits.size(); ++j) {
      RLWE_RETURN_IF_ERROR(curr_a.FusedMulAddInPlace(prev_sub_a_digits[j],
                                                     gk_pads_baby_[j], rns_moduli_));
    }
    ct_pads_baby_.push_back(std::move(curr_a));
    ct_sub_pad_baby_digits_.push_back(std::move(prev_sub_a_digits));
  }

  // std::cout << "Computing diagonals * baby step rotations" << std::endl;
  // Preprocess the databases using the "a" part of Enc(s << i).
  for (auto const& database : databases_) {
    RLWE_RETURN_IF_ERROR(database->PreprocessBSGS(ct_pads_baby_));
  }

  // std::cout << "Last part of preprocessing" << std::endl;
  
  // TODO: Allocate ct_pads_giant_ and giant_digits_ (be careful about the size)
  pad_inner_products_giant_.clear();
  pad_inner_products_giant_.resize(databases_.size());
  ct_pads_giant_.resize(databases_.size());
  ct_sub_pad_giant_digits_.resize(databases_.size());
  int db_ctr = 0;
  for (auto const& database : databases_) {
    std::vector<std::vector<RnsPolynomial>>& inner_product_pads_baby = database->GetBabyInnerProductsPad();

    ct_pads_giant_[db_ctr].resize(inner_product_pads_baby.size());
    ct_sub_pad_giant_digits_[db_ctr].resize(inner_product_pads_baby.size());
    std::vector<RnsPolynomial> ct_pads_final;
    ct_pads_final.reserve(inner_product_pads_baby.size());
    
    for (int i = 0; i < inner_product_pads_baby.size(); i++){
      //std::cout << "Final step of preprocessing " << i << std::endl;
      RLWE_ASSIGN_OR_RETURN(RnsPolynomial pad_acc, RnsPolynomial::CreateZero(log_n, rns_moduli_, true));
      if (inner_product_pads_baby[i].size() != giant_steps_total) {
        std::cout << "=====> ERROR: The size of inner_product_pads_baby[i] not equal to the number of giant steps" << std::endl;
        assert(false);
      }
      //ct_pads_giant_[i].reserve(giant_steps_total);
      //ct_sub_pad_giant_digits_[i].reserve(giant_steps_total);

      // Applying Horner's rule [https://eprint.iacr.org/2018/244.pdf]
      RLWE_RETURN_IF_ERROR(pad_acc.AddInPlace(inner_product_pads_baby[i][giant_steps_total - 1], rns_moduli_));

      for (int Gid = giant_steps_total - 2; Gid >= 0; Gid--){
        //std::cout << "Gid: " << Gid << std::endl;
        // rotate by Gid * baby_steps_total and add to accumulator
        //int rotation_power = small_modexp_2n(5, Gid * baby_steps_total, 1 << log_n);
        //std::cout << "Rotation power: " << rotation_power << std::endl;

        // Rotate pad_acc
        RLWE_ASSIGN_OR_RETURN(RnsPolynomial sub_a,
                              pad_acc.Substitute(giant_step_power, rns_moduli_));
        //std::cout << "Rotation done" << std::endl;

        // Compute gadget inverse, the convert back to NTT form
        if (sub_a.IsNttForm()) {
          RLWE_RETURN_IF_ERROR(sub_a.ConvertToCoeffForm(rns_moduli_));
        }
        RLWE_ASSIGN_OR_RETURN(auto sub_a_digits,
                              rns_gadget_.Decompose(sub_a, rns_moduli_));
        for (auto& digit : sub_a_digits) {
          RLWE_RETURN_IF_ERROR(digit.ConvertToNttForm(rns_moduli_));
        }
        //std::cout << "Gadget decomposition done" << std::endl;

        // Apply galois key
        RLWE_ASSIGN_OR_RETURN(
            pad_acc,
            RnsPolynomial::CreateZero(log_n, rns_moduli_, true));
        for (int j = 0; j < sub_a_digits.size(); ++j) {
          RLWE_RETURN_IF_ERROR(pad_acc.FusedMulAddInPlace(sub_a_digits[j],
                                                        gk_pads_giant_[j], rns_moduli_));
        }
        //std::cout << "Galois key applied" << std::endl;

        ct_pads_giant_[db_ctr][i].push_back(pad_acc);
        ct_sub_pad_giant_digits_[db_ctr][i].push_back(std::move(sub_a_digits));

        // Add inner_product_pads_baby[i][Gid] to pad_acc
        RLWE_RETURN_IF_ERROR(pad_acc.AddInPlace(inner_product_pads_baby[i][Gid], rns_moduli_));


      }
      // Store the accumulator in database object
      ct_pads_final.push_back(std::move(pad_acc));
    }
    
    // Pass the final pads to the database object
    //RLWE_RETURN_IF_ERROR(database->SetPreprocessedFinalPadsBSGS(ct_pads_final));
    
    
    for (int j = 0; j < ct_pads_final.size(); j++) {
      pad_inner_products_giant_[db_ctr].push_back(ct_pads_final[j]);
    }
    db_ctr++;
  }
  
  // std::cout << "Preprocessing done" << std::endl;

  return absl::OkStatus();
}

template <typename RlweInteger>
absl::StatusOr<LinPirResponse> Server<RlweInteger>::HandleRequest(
    const RnsCiphertext& ct_query, const RnsGaloisKey& gk) const {
  // Compute all rotations of the query vector.
  int num_rotations = params_.rows_per_block / 2;
  std::vector<RnsCiphertext> ct_rotated_queries;
  ct_rotated_queries.reserve(num_rotations);
  ct_rotated_queries.push_back(std::move(ct_query));
  for (int i = 1; i < num_rotations; ++i) {
    RLWE_ASSIGN_OR_RETURN(RnsCiphertext ct_sub,
                          ct_rotated_queries[i - 1].Substitute(5));
    RLWE_ASSIGN_OR_RETURN(RnsCiphertext ct_rot, gk.ApplyTo(ct_sub));
    ct_rotated_queries.push_back(std::move(ct_rot));
  }
  // Compute inner products with the databases and serialize.
  LinPirResponse response;
  for (auto const& database : databases_) {
    LinPirResponse::EncryptedInnerProduct inner_product;
    RLWE_ASSIGN_OR_RETURN(std::vector<RnsCiphertext> ct_blocks,
                          database->InnerProductWith(ct_rotated_queries));
    for (auto const& ct : ct_blocks) {
      RLWE_ASSIGN_OR_RETURN(*inner_product.add_ct_blocks(), ct.Serialize());
    }
    *response.add_ct_inner_products() = std::move(inner_product);
  }
  return response;
}

template <typename RlweInteger>
absl::StatusOr<LinPirResponse> Server<RlweInteger>::HandleRequest(
    const ::rlwe::SerializedRnsPolynomial& proto_ct_query_b,
    const google::protobuf::RepeatedPtrField<::rlwe::SerializedRnsPolynomial>&
        proto_gk_key_bs) const {
  // Deserialize the "b" components from request and build the query ciphertext
  // and the Galois key.
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial ct_query_b,
      RnsPolynomial::Deserialize(proto_ct_query_b, rns_moduli_));
  RnsCiphertext ct_query({std::move(ct_query_b), ct_pads_[0]}, rns_moduli_,
                         /*power_of_s=*/1, /*error=*/0, &rns_error_params_,
                         rns_context_);

  std::vector<RnsPolynomial> gk_key_bs;
  gk_key_bs.reserve(proto_gk_key_bs.size());
  for (int i = 0; i < proto_gk_key_bs.size(); ++i) {
    RLWE_ASSIGN_OR_RETURN(
        RnsPolynomial gk_key_b,
        RnsPolynomial::Deserialize(proto_gk_key_bs[i], rns_moduli_));
    gk_key_bs.push_back(std::move(gk_key_b));
  }
  RLWE_ASSIGN_OR_RETURN(
      RnsGaloisKey gk,
      RnsGaloisKey::CreateFromKeyComponents(
          gk_pads_, std::move(gk_key_bs), /*power=*/5, &rns_gadget_,
          rns_moduli_, prng_seed_gk_pad_, params_.prng_type));

  // Compute all rotations of the query vector.
  int num_rotations = params_.rows_per_block / 2;
  std::vector<RnsCiphertext> ct_rotated_queries;
  ct_rotated_queries.reserve(num_rotations);
  ct_rotated_queries.push_back(std::move(ct_query));
  // TODO TODO CAUTION: since giant_factor * baby_factor is > #rotations, be careful because the #rotations will get exhausted earlier than expected
  for (int i = 1; i < num_rotations; ++i) {
    RLWE_ASSIGN_OR_RETURN(RnsCiphertext ct_sub,
                          ct_rotated_queries[i - 1].Substitute(5));
    RLWE_ASSIGN_OR_RETURN(RnsCiphertext ct_rot,
                          gk.ApplyToWithRandomPad(
                              ct_sub, ct_sub_pad_digits_[i - 1], ct_pads_[i]));
    ct_rotated_queries.push_back(std::move(ct_rot));
  }

  // Compute inner products with the databases and serialize them.
  LinPirResponse response;
  response.mutable_ct_inner_products()->Reserve(databases_.size());
  for (auto const& database : databases_) {
    RLWE_ASSIGN_OR_RETURN(
        std::vector<RnsCiphertext> ct_blocks,
        database->InnerProductWithPreprocessedPads(ct_rotated_queries));
    LinPirResponse::EncryptedInnerProduct inner_product;
    inner_product.mutable_ct_blocks()->Reserve(ct_blocks.size());
    for (auto const& ct : ct_blocks) {
      RLWE_ASSIGN_OR_RETURN(*inner_product.add_ct_blocks(), ct.Serialize());
    }
    *response.add_ct_inner_products() = std::move(inner_product);
  }
  return response;
}

template <typename RlweInteger>
absl::StatusOr<LinPirResponse> Server<RlweInteger>::HandleRequestBSGS(
    const ::rlwe::SerializedRnsPolynomial& proto_ct_query_b,
    const google::protobuf::RepeatedPtrField<::rlwe::SerializedRnsPolynomial>&
        proto_gk_key_bs_baby,
    const google::protobuf::RepeatedPtrField<::rlwe::SerializedRnsPolynomial>&
        proto_gk_key_bs_giant) const {
  
  // Deserialize the "b" components from request and build the query ciphertext
  // and the Galois key.
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial ct_query_b,
      RnsPolynomial::Deserialize(proto_ct_query_b, rns_moduli_));
  RnsCiphertext ct_query({std::move(ct_query_b), ct_pads_baby_[0]}, rns_moduli_,
                         /*power_of_s=*/1, /*error=*/0, &rns_error_params_,
                         rns_context_);

  //std::cout << "Deserialized query polynomial" << std::endl;
  int num_slots = 1 << params_.log_n;
  std::vector<RlweInteger> one_vector(num_slots, 0);
  RLWE_ASSIGN_OR_RETURN(auto 
  encoder, rlwe::FiniteFieldEncoder<ModularInt>::Create(rns_context_));
  RLWE_ASSIGN_OR_RETURN(
            RnsPolynomial one,
            encoder.EncodeBfv(one_vector, rns_moduli_, /*is_scaled=*/false));

  std::vector<RnsPolynomial> gk_key_bs_baby, gk_key_bs_giant;
  gk_key_bs_baby.reserve(proto_gk_key_bs_baby.size());
  gk_key_bs_giant.reserve(proto_gk_key_bs_giant.size());
  for (int i = 0; i < proto_gk_key_bs_baby.size(); ++i) {
    RLWE_ASSIGN_OR_RETURN(
        RnsPolynomial gk_key_b,
        RnsPolynomial::Deserialize(proto_gk_key_bs_baby[i], rns_moduli_));
    gk_key_bs_baby.push_back(std::move(gk_key_b));
  }
  for (int i = 0; i < proto_gk_key_bs_giant.size(); ++i) {
    RLWE_ASSIGN_OR_RETURN(
        RnsPolynomial gk_key_b,
        RnsPolynomial::Deserialize(proto_gk_key_bs_giant[i], rns_moduli_));
    gk_key_bs_giant.push_back(std::move(gk_key_b));
  }
  
  int num_blocks = DivideAndRoundUp(params_.db_rows, params_.rows_per_block);
  int num_rotations = params_.rows_per_block / 2;
  auto baby_giant = find_optimal_bsgs(num_rotations, num_blocks);
  int baby_steps_total = baby_giant.first;
  int giant_steps_total = baby_giant.second;
  //int giant_steps_total = std::ceil(std::sqrt(num_rotations) / params_.baby_steps_bias_factor);
  //int baby_steps_total = std::ceil((float)num_rotations / giant_steps_total);

  int baby_step_power = 5;
  int giant_step_power = small_modexp_2n(baby_step_power, baby_steps_total, 1 << params_.log_n); 

  RLWE_ASSIGN_OR_RETURN(
      RnsGaloisKey gk_baby,
      RnsGaloisKey::CreateFromKeyComponents(
          gk_pads_baby_, std::move(gk_key_bs_baby), /*power=*/baby_step_power, &rns_gadget_,
          rns_moduli_, prng_seed_gk_pad_baby_, params_.prng_type));

  RLWE_ASSIGN_OR_RETURN(
      RnsGaloisKey gk_giant,
      RnsGaloisKey::CreateFromKeyComponents(
          gk_pads_giant_, std::move(gk_key_bs_giant), /*power=*/giant_step_power, &rns_gadget_,
          rns_moduli_, prng_seed_gk_pad_giant_, params_.prng_type));

  // Compute all rotations of the query vector.
  //std::cout << "Computing baby step rotations" << std::endl;
  std::vector<RnsCiphertext> ct_rotated_queries_baby;
  ct_rotated_queries_baby.reserve(baby_steps_total);
  ct_rotated_queries_baby.push_back(std::move(ct_query));

  for (int i = 1; i < baby_steps_total; ++i) {
    RLWE_ASSIGN_OR_RETURN(RnsCiphertext ct_sub,
                          ct_rotated_queries_baby[i - 1].Substitute(5));
    RLWE_ASSIGN_OR_RETURN(RnsCiphertext ct_rot,
                          gk_baby.ApplyToWithRandomPad(
                              ct_sub, ct_sub_pad_baby_digits_[i - 1], ct_pads_baby_[i]));
    ct_rotated_queries_baby.push_back(std::move(ct_rot));
  }

  // Compute inner products with the databases and serialize them.
  LinPirResponse response;
  response.mutable_ct_inner_products()->Reserve(databases_.size());
  int db_ctr = 0;
  
  for (auto const& database : databases_) {
    std::vector<RnsCiphertext> ct_blocks;

    RLWE_ASSIGN_OR_RETURN(
        std::vector<std::vector<RnsCiphertext>> ct_inner_prods,
        database->InnerProductWithPreprocessedPadsBSGS(ct_rotated_queries_baby));
    //std::cout << "Computed inner products (1st stage of BSGS)" << std::endl;

    for (int i = 0; i < ct_inner_prods.size(); i++) {
      if (ct_inner_prods[i].size() != giant_steps_total) {
        std::cout << "=====> ERROR: The size of ct_inner_prods[i] not equal to the number of giant steps" << std::endl;
        assert(false);
      }
      // Horner's method [https://eprint.iacr.org/2018/244.pdf]

      //auto error_params = ct_inner_prods[0][0].ErrorParams();
      //RnsCiphertext ct_acc(
      //  RnsCiphertext::CreateZero(rns_moduli_, error_params));
      
      //RLWE_RETURN_IF_ERROR(
      //      ct_acc.FusedAbsorbAddInPlace(
      //          ct_inner_prods[i][giant_steps_total - 1], one));
      //std::cout << "Adding last (which is the first) inner product" << std::endl;
      //RLWE_RETURN_IF_ERROR(
      //      ct_acc.AddInPlace(
      //          ct_inner_prods[i][giant_steps_total - 1]));
      
      auto ct_acc = ct_inner_prods[i][giant_steps_total - 1];
      //if (db_ctr > 0) {
      //  ct_acc = ct_inner_prods[i][giant_steps_total - 1];
      //}
      int ctr = 0;
      for (int Gid = giant_steps_total - 2; Gid >= 0; Gid--) {
        //std::cout << "Gid: " << Gid << std::endl;
        //int rotation_power = small_modexp_2n(5, Gid * baby_steps_total, 1 << log_n);

        RLWE_ASSIGN_OR_RETURN(RnsCiphertext ct_sub,
                              ct_acc.Substitute(giant_step_power));
        RLWE_ASSIGN_OR_RETURN(ct_acc,
                              gk_giant.ApplyToWithRandomPad(
                                  ct_sub, ct_sub_pad_giant_digits_[db_ctr][i][ctr], ct_pads_giant_[db_ctr][i][ctr]));
        ctr++;
        RLWE_RETURN_IF_ERROR(ct_acc.AddInPlace(ct_inner_prods[i][Gid]));
        //RLWE_RETURN_IF_ERROR(
        //    ct_acc.FusedAbsorbAddInPlace(
        //        ct_inner_prods[i][Gid], one));

      }
      RLWE_RETURN_IF_ERROR(ct_acc.SetPadComponent(pad_inner_products_giant_[db_ctr][i]));
      ct_blocks.push_back(std::move(ct_acc));
    }


    LinPirResponse::EncryptedInnerProduct inner_product;
    inner_product.mutable_ct_blocks()->Reserve(ct_blocks.size());
    for (auto const& ct : ct_blocks) {
      RLWE_ASSIGN_OR_RETURN(*inner_product.add_ct_blocks(), ct.Serialize());
    }

    *response.add_ct_inner_products() = std::move(inner_product);
    db_ctr++;
  }
  return response;
}

template class Server<Uint32>;
template class Server<Uint64>;

}  // namespace linpir
}  // namespace hintless_pir
