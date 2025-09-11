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

#include "hintless_simplepir/simplepir.h"
#include "verisimplepir/src/lib/pir/preproc_pir.h"
#include "verisimplepir/src/lib/pir/pir.h"

#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <iostream>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "lwe/types.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "lwe/matrix_conversion.h"
#include "hintless_simplepir/client.h"
#include "hintless_simplepir/database_hwy.h"
#include "hintless_simplepir/parameters.h"
#include "hintless_simplepir/server.h"
#include "linpir/parameters.h"

namespace hintless_pir {

namespace hintless_simplepir {
namespace {

using RlweInteger = Parameters::RlweInteger;

const int rows_db = 2048;
const int cols_db = 1024;
bool use_static_db = true;

const Parameters kParameters{
    .db_rows = rows_db,
    .db_cols = cols_db,
    .db_record_bit_size = 8, // Don't include stacking here.
    .lwe_secret_dim = 1408,
    .offline_lwe_secret_dim = 2816,
    .lwe_modulus_bit_size = 32,
    .lwe_plaintext_bit_size = 8, // <= 8
    .db_stack_cells = 1,
    .lwe_error_variance = 8,
    .linpir_params =
        linpir::RlweParameters<RlweInteger>{
            .log_n = 12, 
            .qs = {72057594021150721ULL, 72057589742960641ULL}, 
#ifdef SH_RUN
            .ts = {2056193, 1990657},
#else //SH_RUN 
            .ts = {SOL_NTT_PRIME},
#endif //SH_RUN
            .gadget_log_bs = {16, 16}, 
            .error_variance = 8,
            .prng_type = rlwe::PRNG_TYPE_HKDF,
#ifdef BSGS
            .rows_per_block = std::min(rows_db, 2048),
#else //BSGS
            .rows_per_block = 1024, // From HintlessPIR
#endif //BSGS
            .db_rows = rows_db, // Need for BSGS optimality computation
        },
    .prng_type = rlwe::PRNG_TYPE_HKDF,
};

// Client's query index. This is database entry that the client will access via PIR
uint64_t index = cols_db - 4; // 4-th column from the right starting at 0th row (indexing from 0).

// All recorded metrics
#define HINT_MB "Hints (MiB)"
#define ONLINE_STATE_KB "Online State (KiB)"
#define OFFLINE_UP_KB "Offline Up (KiB)"
#define OFFLINE_DOWN_KB "Offline Down (KiB)"
#define PREPARE_UP_KB "Prep Up (Kib)"
#define PREPARE_DOWN_KB "Prep Down (KiB)"
#define QUERY_UP_KB "Query Up (KiB)"
#define QUERY_DOWN_KB "Query Down (KiB)"
#define GLOBAL_PREPR_S "Global Prepr (s)"
#define SV_CLIENT_SPEC_PREPR_S "Server Per-Client Prepr (s)"
#define CLIENT_PREPR_S "Client Local Prepr (s)"
#define CLIENT_PREP_PRE_S "Client Prepa Pre Req (s)"
#define SERVER_PREP_S "Server Prepa Comp (s)"
#define CLIENT_PREP_POST_S "Client Prepa Post Req (s)"
#define CLIENT_Q_REQ_GEN_MS "Query: Client Req Gen (ms)"
#define SERVER_Q_RESP_S "Query: Server Comp (s)"
#define CLIENT_Q_DEC_MS "Query: Client Decryption (ms)"
#define CLIENT_Q_VER_MS "Query: Client Verification (ms)"
#define LONG_TERM_STATE_KB "Long Term State (KiB)"

TEST(HintlessSimplePir, EndToEndPIRTest) {
  // Some calculation around client's query index given how database is laid out
  uint64_t row_idx_start = index / cols_db;
  if (row_idx_start % kParameters.db_stack_cells != 0) {
    std::cout << "\n\nError: row index not divisible by stack cells\n\n" << std::endl;
    assert(false);
  }
  uint64_t row_idx_end = row_idx_start + kParameters.db_stack_cells;
  uint64_t col_idx = index % cols_db;

  // Important metrics
  map<string, double> dict;
  dict[HINT_MB] = 0; // mb is MiB, kb is KiB, s is sec, ms is ms
  dict[ONLINE_STATE_KB] = 0;
  dict[OFFLINE_UP_KB] = 0;
  dict[OFFLINE_DOWN_KB] = 0;
  dict[PREPARE_UP_KB] = 0;
  dict[PREPARE_DOWN_KB] = 0;
  dict[QUERY_UP_KB] = 0;
  dict[QUERY_DOWN_KB] = 0;
  dict[GLOBAL_PREPR_S] = 0;
  dict[SV_CLIENT_SPEC_PREPR_S] = 0;
  dict[CLIENT_PREPR_S] = 0;
  dict[CLIENT_PREP_PRE_S] = 0;
  dict[SERVER_PREP_S] = 0;
  dict[CLIENT_PREP_POST_S] = 0;
  dict[CLIENT_Q_REQ_GEN_MS] = 0;
  dict[SERVER_Q_RESP_S] = 0;
  dict[CLIENT_Q_DEC_MS] = 0;
  dict[CLIENT_Q_VER_MS] = 0;
  dict[LONG_TERM_STATE_KB] = 0;
  
  const uint64_t row_in = kParameters.db_rows;
  const uint64_t col_in = kParameters.db_cols;
  const uint64_t lattice_n_in = kParameters.lwe_secret_dim;
  const uint64_t offline_lattice_n_in = kParameters.offline_lwe_secret_dim;
  const uint64_t d_in = kParameters.db_record_bit_size;
  const uint64_t N = row_in * col_in;
  const uint64_t logp_in = kParameters.lwe_plaintext_bit_size;
  
  std::cout << "Input params: N = " << N << " d = " << d_in << std::endl;
  std::cout << "database size: " << N*d_in / (8.0*(1ULL<<30)) << " GiB\n";
  std::cout << "Rows per block: " << kParameters.linpir_params.rows_per_block << std::endl;
  std::cout << "Shards: " << kParameters.db_record_bit_size / kParameters.lwe_plaintext_bit_size << " and Stacks: " << kParameters.db_stack_cells << std::endl;

#ifdef BSGS
  std::cout << "-----   Using BSGS   -----\n" << std::endl;
  if (kParameters.db_rows % kParameters.linpir_params.rows_per_block != 0) {
    std::cout << "ERROR: rows_per_block must divide db_rows" << std::endl;
    assert(false);
  }
#endif //BSGS

#ifdef FAKE_RUN
  std::cout << "-----   Fake Run   -----" << std::endl;
#endif //FAKE_RUN

#ifdef SH_RUN
  std::cout << "-----   Semi-Honest Run   -----" << std::endl;
#endif //SH_RUN

  double start, start2, start_prep, start_on, end, end2, end_prep, end_on;

	// Create VSPIR Object
  VeriSimplePIR pir(row_in, col_in, lattice_n_in, d_in, logp_in, true, false, false, false, 1, true, false);
  //std::cout << "Database params: "; pir.dbParams.print();

  // Create server and fill in random database records.
  std::unique_ptr<Server> server;

  start = currentDateTime();
#ifdef FAKE_RUN
  std::cout << "Using fake DB; correctness won't hold; this is only for evaluation" << std::endl;
  ASSERT_OK_AND_ASSIGN(server,
                       Server::CreateWithRandomDatabaseRecordsFake(kParameters));
#else //FAKE_RUN
  if (use_static_db) {
    std::cout << "Using static DB" << std::endl;
    ASSERT_OK_AND_ASSIGN(server,
                       Server::CreateWithFastDatabaseRecords(kParameters));
  }
  else { // use random DB
    std::cout << "Using random DB" << std::endl;
    ASSERT_OK_AND_ASSIGN(server,
                       Server::CreateWithRandomDatabaseRecords(kParameters));
  }
#endif //FAKE_RUN
  end = currentDateTime();
  //std::cout << "[==> TIMER  <==] Server database creation and transpose time: " << (end-start) << " ms | " << (end-start)/1000 << " sec" << std::endl;
  
  const Database* database = server->GetDatabase();
  // Preprocess the server and get public parameters.
  start_prep = currentDateTime();
  start = currentDateTime();
#ifdef BSGS
  ASSERT_OK(server->PreprocessBSGS());
#else //BSGS
  ASSERT_OK(server->Preprocess());
#endif //BSGS
  end = currentDateTime();
  // Server global preprocessing #1 (only main hint, i.e. overall LWE instance, and RLWE part) time
  dict[GLOBAL_PREPR_S] += (end-start)/1000;

  auto public_params = server->GetPublicParams();

  // Create a client and issue request.
  start = currentDateTime();
  ASSERT_OK_AND_ASSIGN(auto client, Client::Create(kParameters, public_params));
  end = currentDateTime();
  //std::cout << RED << "[==> TIMER  <==] Client creation time: " << (end-start) << " ms | " << (end-start)/1000 << " sec" << END << std::endl;

  // This code moves back and forth between HintlessPIR code and VeriSimplePIR (VSPIR) code. So some type-casting is needed
  // Get A, H for VSPIR

  start = currentDateTime();
  // Get D, A, H from Hintless
  int num_shards = database->Data().size();
  auto D_raw_vec = database->Data();
  std::vector<lwe::Matrix> D_vec;
  D_vec.reserve(num_shards);
  int bits_for_db_mask = kParameters.db_record_bit_size;
  if (bits_for_db_mask > kParameters.lwe_plaintext_bit_size) {
    bits_for_db_mask = kParameters.lwe_plaintext_bit_size; // smaller of the two
  }
  for (int i = 0; i < num_shards; i++) {
    lwe::Matrix D = ExportRawMatrix(D_raw_vec[i], kParameters.db_rows, bits_for_db_mask);
    D_vec.push_back(D);
  }
  const lwe::Matrix* A1_ref = database->LweQueryPad();
  auto A1 = ImportLweMatrix(*A1_ref);
  auto H_vec = database->Hints();

  // Convert D, A, H to VSPIR //
	// A
	VMatrix A_vsp(col_in, kParameters.lwe_secret_dim);
  convertToVeriSimplePIRMatrix(A_vsp, A1);

	// H (directly rather than recomputing the hint)
  std::vector<VMatrix> H_vsp_vec;
  H_vsp_vec.reserve(num_shards);
  for (int i = 0; i < num_shards; i++) {
    VMatrix H_vsp(H_vec[0].size(), H_vec[0][0].size());
	  convertToVeriSimplePIRMatrixRowFormIn(H_vsp, H_vec[i]);
    H_vsp_vec.push_back(H_vsp);
  }
  // Main (overall LWE instance) hint size (MB)
  dict[HINT_MB] += (H_vsp_vec.size() * H_vsp_vec[0].rows * H_vsp_vec[0].cols * sizeof(Elem))/(1ULL << 20);
  end = currentDateTime();
  // Server converting Hintless -> VeriSimplePIR time
  dict[GLOBAL_PREPR_S] += (end-start)/1000;
  
	// Conversions done

#ifndef SH_RUN // Only run proof preprocessing if SH_RUN is not set

  const BinaryMatrix C = pir.PreprocSampleC();

  start = currentDateTime();
  const Matrix64 A_2 = pir.PreprocInitNew64();
	auto A_2_hintless = convertToHintlessPIRMatrixColumnForm64(A_2, row_in, offline_lattice_n_in);
  end = currentDateTime();
  // Expanding A_2 (A for offline-only LWE instance) time
  dict[GLOBAL_PREPR_S] += (end-start)/1000;
  dict[CLIENT_PREPR_S] += (end-start)/1000;
  
  // To get time to expand main A for client's preprocessing time. Note that A_2 is for offline LWE instance which is used to compute Z = C*D only.
  start = currentDateTime();
  const VMatrix A_temp = pir.PreprocInitNew(lattice_n_in);
  auto A_temp_hintless = convertToHintlessPIRMatrixColumnForm(A_temp, col_in, lattice_n_in);
  end = currentDateTime();
  // Expanding main A time
  dict[CLIENT_PREPR_S] += (end-start)/1000; // Already included in server time from above

  //std::cout << "Computing H_2 = D^T * A_2" << std::endl;
  start = currentDateTime();
#ifdef FAKE_RUN
  auto H_2_hintless_vec = server->TransposeProductFake64();
#else //FAKE_RUN
	ASSERT_OK_AND_ASSIGN(auto H_2_hintless_vec, server->TransposeProduct64(A_2_hintless));
#endif //FAKE_RUN

  std::vector<Matrix64> H_2_vec;
  H_2_vec.reserve(num_shards);
  for (int i = 0; i < num_shards; i++) {
    Matrix64 H_2(col_in, offline_lattice_n_in);
	  convertToVeriSimplePIRMatrixRowFormIn64(H_2, H_2_hintless_vec[i]);
    H_2_vec.push_back(H_2);
  }
  end = currentDateTime();
  // Server global H_2 (hint for offline-only LWE instance) compute time
  dict[GLOBAL_PREPR_S] += (end-start)/1000;
  // Large offline-only LWE's hint H_2 size (MB)
  dict[HINT_MB] += (H_2_vec.size() * H_2_vec[0].rows * H_2_vec[0].cols * sizeof(Elem64))/(1ULL << 20);

  //std::cout << "Encrypting challenge C for preprocessing phase" << std::endl;
  start = currentDateTime();
  const auto preproc_ct_sk_pair = pir.PreprocClientMessageNew64(A_2, C);

  const auto preproc_cts = std::get<0>(preproc_ct_sk_pair);
  const auto preproc_sks = std::get<1>(preproc_ct_sk_pair);
  end = currentDateTime();
  // Client encrypting C time
  dict[CLIENT_PREPR_S] += (end-start)/1000;
  // Client offline upload size (KB)
  dict[OFFLINE_UP_KB] += (preproc_cts.size() * preproc_cts[0].rows * preproc_cts[0].cols * sizeof(Elem64))/(1ULL << 10);

  //std::cout << "Computing v_i = D^T * u_i, for all i, in preprocessing" << std::endl;
  start = currentDateTime();
	auto preproc_cts_hintless = convertToHintlessPIRMatricesFromVectors64(preproc_cts, row_in);
	ASSERT_OK_AND_ASSIGN(auto D_T_ct, server->TransposeProductCiphertexts64(preproc_cts_hintless));
	
  std::vector<std::vector<Matrix64>> preproc_res_cts_vec;
  preproc_res_cts_vec.reserve(num_shards);
  for (int i = 0; i < num_shards; i++) {
    auto preproc_res_cts = convertToVSPIRPIRVectors64(D_T_ct[i]);
    preproc_res_cts_vec.push_back(preproc_res_cts);
  }
  end = currentDateTime();
  // Server (client-specific) D^T*C^T time
  dict[SV_CLIENT_SPEC_PREPR_S] += (end-start)/1000;
  // Server offline response size (KB)
  dict[OFFLINE_DOWN_KB] += (preproc_res_cts_vec.size() * preproc_res_cts_vec[0].size() * preproc_res_cts_vec[0][0].rows * preproc_res_cts_vec[0][0].cols * sizeof(Elem64))/(1ULL << 10);

	//std::cout << "Decrypting Z in preprocessing phase" << std::endl;
  start = currentDateTime();
  std::vector<VMatrix> Z_vec;
  Z_vec.reserve(num_shards);
  for (int i = 0; i < num_shards; i++) {
    VMatrix Z = pir.PreprocRecoverZNew64(H_2_vec[i], preproc_sks, preproc_res_cts_vec[i], row_in);
    Z_vec.push_back(Z);
  }
  end = currentDateTime();
  // Client decrypting Z time
  dict[CLIENT_PREPR_S] += (end-start)/1000;

  //std::cout << "Sampling random vector for proof compression" << std::endl;
  start = currentDateTime();
  int rlc_count = DivideAndRoundUp(STAT_SEC_PARAM, LOG_Q);
  std::cout << "Proof compression: no. of rows in Z_prime = RC = " << rlc_count << std::endl;
  VMatrix R(rlc_count, STAT_SEC_PARAM);
  random(R, SOL_NTT_PRIME);

  //std::cout << "Compressing the challenge C" << std::endl;
  Matrix RC = matMulRightBinary(R, C);

  //std::cout << "Compressing the proof Z" << std::endl;
  std::vector<VMatrix> RZ_vec;
  RZ_vec.reserve(num_shards);
  for (int i = 0; i < num_shards; i++) {
    VMatrix RZ = matMul(R, Z_vec[i]);
    RZ_vec.push_back(RZ);
  }
  end = currentDateTime();
  // Client compression of Z time
  dict[CLIENT_PREPR_S] += (end-start)/1000;

  // Client compressed ZA=CH check
  start = currentDateTime();
  for (int i = 0; i < num_shards; i++) {
    //std::cout << "(Compressed) verifying shard " << i << " with primary A, H" << std::endl;
    pir.VerifyPreprocZCompressed(RZ_vec[i], A_vsp, RC, H_vsp_vec[i]);
  }
  end = currentDateTime();
  // Client offline verification time
  dict[CLIENT_PREPR_S] += (end-start)/1000;
  std::cout << "Success: (Compressed) Z in preprocessing verified and preprocessing phase done.\n";
  
  end_prep = currentDateTime();
  std::cout << "[==> TIMER  <==] Total preprocessing (server + client compute) time: " << (end_prep-start_prep) << " ms | " << (end_prep-start_prep)/1000 << " sec" << std::endl;

#endif //SH_RUN

  // Long-term (persistent) state only includes the compressed proof RZ. On the other hand, RC is expandable when needed.
  dict[LONG_TERM_STATE_KB] += (RZ_vec.size() * RZ_vec[0].rows * RZ_vec[0].cols * sizeof(Elem))/(1ULL << 10);

  // Query-level preprocessing/prepare + online phase

  start_on = currentDateTime(); // Captures time of query-level prepr + online

  // Prepare phase = query-level preprocessing
  double start_prepare, end_prepare;

  start_prepare = currentDateTime();
	// Generate HintlessPIR client request + output sk from Hintless to VSPIR
  start = currentDateTime();
  ASSERT_OK_AND_ASSIGN(auto As_s_tuple, client->Compute_A_times_s());
  end = currentDateTime();
  // Client expansion of A and A*s computation time
  dict[CLIENT_PREP_PRE_S] += (end-start)/1000;

  auto As_vsp = std::get<0>(As_s_tuple);
  auto As = std::get<1>(As_s_tuple);
  auto s_lwe = std::get<2>(As_s_tuple);

  // Size of of A*s
  dict[ONLINE_STATE_KB] += (As_vsp.rows * As_vsp.cols * sizeof(Elem))/(1ULL << 10);
  // Size of of s_lwe
  dict[ONLINE_STATE_KB] += (s_lwe.Key().size() * sizeof(Elem))/(1ULL << 10); 

  // Prepare RLWE parts
  // Client's RLWE request
  start = currentDateTime();
#ifdef BSGS
  ASSERT_OK_AND_ASSIGN(auto prepare_req, client->PrepareLinPirGivenS_BSGS(s_lwe));
#else //BSGS
  ASSERT_OK_AND_ASSIGN(auto prepare_req, client->PrepareLinPirGivenS(s_lwe));
#endif //BSGS
  end = currentDateTime();
  // Client prepare req gen time
  dict[CLIENT_PREP_PRE_S] += (end-start)/1000;
  // Client prepare request KB
  dict[PREPARE_UP_KB] += (prepare_req.ByteSizeLong() / 1024);
  // Client prepare request (rotation keys, RLWE cipher of s_lwe) size
  dict[ONLINE_STATE_KB] += (prepare_req.ByteSizeLong() / 1024);

  // Server's RLWE response
  start = currentDateTime();
#ifdef BSGS
  ASSERT_OK_AND_ASSIGN(auto prepare_response, server->HandlePrepareRequestBSGS(prepare_req)); 
#else //BSGS
  ASSERT_OK_AND_ASSIGN(auto prepare_response, server->HandlePrepareRequest(prepare_req));
#endif //BSGS
  end = currentDateTime();
  // Server prepare response time (LinPir H*s compute time)
  dict[SERVER_PREP_S] += (end-start)/1000;
  // Server prepare response KB
  dict[PREPARE_DOWN_KB] += (prepare_response.ByteSizeLong() / 1024);
  // Client prepare response from server (Hs) size
  dict[ONLINE_STATE_KB] += (prepare_response.ByteSizeLong() / 1024);

  // Client decrypt RLWE response
  start = currentDateTime();
  ASSERT_OK_AND_ASSIGN(auto hs_vec, client->RecoverHsPreparePhase(prepare_response));
  end = currentDateTime();
  // Client prepare recover time (LinPir response decryption)
  dict[CLIENT_PREP_POST_S] += (end-start)/1000;

#ifndef SH_RUN // Only run verification if SH_RUN is not set
  // Size of Z' = RZ
  dict[ONLINE_STATE_KB] += (RZ_vec.size() * RZ_vec[0].rows * RZ_vec[0].cols * sizeof(Elem))/(1ULL << 10);
  // Size of C' = RC 
  dict[ONLINE_STATE_KB] += (RC.rows * RC.cols * sizeof(Elem))/(1ULL << 10);

  //std::cout << "Doing verification: H*s as part of prepare phase" << std::endl;
  start2 = currentDateTime();
  for (int i = 0; i < num_shards; i++) {
    VMatrix h_times_s = convertToVeriSimplePIRColumnMatrixSized(hs_vec[i], row_in);
    // Size of decrypted H*s (i.e. w)
    dict[ONLINE_STATE_KB] += (h_times_s.rows * h_times_s.cols * sizeof(Elem))/(1ULL << 10);
    // Verify
    pir.HsVerifyCompressed(As_vsp, h_times_s, RZ_vec[i], RC);
  }
  end2 = currentDateTime();
  // Client prepare verification of H*s time
  dict[CLIENT_PREP_POST_S] += (end2-start2)/1000;
#endif //SH_RUN

  end_prepare = currentDateTime();
  std::cout << BLUE << "[==> TIMER  <==] Prepare phase (query-level preprocessing) total (client + server) time: " << (end_prepare-start_prepare) << " ms | " << (end_prepare-start_prepare)/1000 << " sec" << END << std::endl;

  double start_online, end_online;
  start_online = currentDateTime();
  lwe::Vector query_vsp = lwe::Vector::Zero(col_in);
  
  // Client online request gen
  start = currentDateTime();
  ASSERT_OK_AND_ASSIGN(auto request_and_query, client->GenerateRequestGivenAsSkipLinPir(index, As, s_lwe));
  end = currentDateTime();
  // Client online request generation time
  dict[CLIENT_Q_REQ_GEN_MS] += (end-start);
	
  auto request = request_and_query.first;
	auto query = request_and_query.second;
  // Client online request (LWE only) KB
  dict[QUERY_UP_KB] += (request.ByteSizeLong() / 1024);
  // Client request (LWE query only) size
  dict[ONLINE_STATE_KB] += (request.ByteSizeLong() / 1024);

	VMatrix query_vsp__ = convertToVeriSimplePIRColumnMatrix(query);

	// Server handles the online request
  start2 = currentDateTime();
  ASSERT_OK_AND_ASSIGN(auto response, server->HandleRequestSkipLinPir(request));
  end2 = currentDateTime();
  // Server-only online time (only D*u)
  dict[SERVER_Q_RESP_S] += (end2-start2)/1000;
  // Server online response KB 
  dict[QUERY_DOWN_KB] += (response.ByteSizeLong() / 1024);
  // Client response from server (D*u) size
  dict[ONLINE_STATE_KB] += (response.ByteSizeLong() / 1024);

  // Client decrypts the online response
  start = currentDateTime();
  ASSERT_OK_AND_ASSIGN(auto rec_v, client->RecoverRecordAndGiveV(response, hs_vec));
  end = currentDateTime();
  // Client record recovery time (given precomputed Hs)
  dict[CLIENT_Q_DEC_MS] += (end-start);

  auto record = std::get<0>(rec_v);
	auto resp_vec = std::get<1>(rec_v);

#ifndef SH_RUN // Only run verification if SH_RUN is not set
  start = currentDateTime();
  //std::cout << "Doing verification: D*u" << std::endl;
  for (int i = 0; i < num_shards; i++) {
    auto response_vsp__ = convertToVeriSimplePIRColumnMatrixStdVec(resp_vec[i]);
    pir.PreVerifyCompressed(query_vsp__, response_vsp__, RZ_vec[i], RC);
  }
  //std::cout << "Online verification done for all shards" << std::endl;
  end = currentDateTime();
  // Client online verification time (only D*u)
  dict[CLIENT_Q_VER_MS] += (end-start);
#endif //SH_RUN

  end_on = currentDateTime();
  std::cout << BLUE << "[==> TIMER  <==] Query-level preprocessing + online time: " << (end_on-start_on) << " ms | " << (end_on-start_on)/1000 << " sec" << END << std::endl;

  end_online = currentDateTime();
  std::cout << BLUE << "[==> TIMER  <==] Online-only (w/o prepare) phase total (client + server) time: " << (end_online-start_online) << " ms | " << (end_online-start_online)/1000 << " sec" << END << std::endl;
  
  if (kParameters.db_stack_cells > 1) { // Stacking
    ASSERT_OK_AND_ASSIGN(auto expected, database->RecordWithStacking(index));
    EXPECT_EQ(record, expected);
  }
  else {
    ASSERT_OK_AND_ASSIGN(auto expected, database->Record(index));
    EXPECT_EQ(record, expected);
  }
  std::cout << "End-to-end test done." << std::endl;
  
  std::cout << "----------------------------------" << std::endl;
  std::cout << "Recorded metrics" << std::endl;
  std::cout << "Database size: " << N*d_in / (8.0*(1ULL<<30)) << " GiB\n";
  std::cout << "Shards: " << kParameters.db_record_bit_size / 
                                 kParameters.lwe_plaintext_bit_size << " and Stacks: " << kParameters.db_stack_cells << std::endl;
  for (auto it : dict) {
    std::cout << it.first << " : " << it.second << std::endl;
  }

}

}  // namespace
}  // namespace hintless_simplepir

}  // namespace hintless_pir
