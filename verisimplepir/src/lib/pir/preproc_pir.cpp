#include "preproc_pir.h"
#include <cmath>
#include <functional>

// Operations for preprocessing phase

// Sample A_2 for the preprocessing
Multi_Limb_Matrix VeriSimplePIR::PreprocInit() const {
    return preproc_lhe.genPublicA(ell);
};

Matrix VeriSimplePIR::PreprocInitNew() const {
    return lhe.genPublicA(ell);
};

Matrix VeriSimplePIR::PreprocInitNew(uint64_t lwe_dim) const {
    return lhe.genPublicA(m, lwe_dim);
};

Matrix64 VeriSimplePIR::PreprocInitNew64() const {
    return lhe.genPublicA64(ell);
};

Matrix VeriSimplePIR::PreprocInitNewSecondary() const {
    return lhe.genPublicASecondary(ell);
};

Matrix VeriSimplePIR::PreprocInitNewSecondary(uint64_t lwe_dim) const {
    return lhe.genPublicASecondary(m, lwe_dim);
};

Multi_Limb_Matrix VeriSimplePIR::PreprocFakeInit() const {
    Multi_Limb_Matrix A(ell, lhe.n);
    random_fast(A.q_data);
    random_fast(A.kappa_data, preproc_lhe.kappa);
    return A;
}

Multi_Limb_Matrix VeriSimplePIR::PreprocGenerateHint(const Multi_Limb_Matrix& A, const Matrix& D) const {
    if (D.rows != m || D.cols != ell) {
        std::cout << "database dimension mismatch! input should be D^T\n";
        assert(false);
    }

    Multi_Limb_Matrix H = matMul(D, A, preproc_lhe.kappa);
    return H;
}

Matrix VeriSimplePIR::PreprocGenerateHintNew(const Matrix& A, const Matrix& D) const {
    if (D.rows != m || D.cols != ell) {
        std::cout << "database dimension mismatch! input should be D^T\n";
        assert(false);
    }

    Matrix H = matMul(D, A);
    return H;
}

Multi_Limb_Matrix VeriSimplePIR::PreprocGenerateFakeHint() const {
    Multi_Limb_Matrix H(m, lhe.n);
    random_fast(H.q_data);
    random_fast(H.kappa_data, preproc_lhe.kappa);
    return H;
}

// samples the plaintext C to be encrypted
BinaryMatrix VeriSimplePIR::PreprocSampleC() const {
    BinaryMatrix C(stat_sec_param, ell);
    random(C);
    // constant(C, 0); std::cout << "change me back\n";
    return C;
}

// samples the plaintext C to be encrypted
BinaryMatrix VeriSimplePIR::PreprocSampleCExplicitDim(uint32_t cols) const {
    BinaryMatrix C(stat_sec_param, cols);
    random(C);
    // constant(C, 0); std::cout << "change me back\n";
    return C;
}

// Encrypts C with a fresh key. output is list of ciphertexts and secret keys pair  
std::pair<std::vector<Multi_Limb_Matrix>, std::vector<Multi_Limb_Matrix>>
VeriSimplePIR::PreprocClientMessage(
    const Multi_Limb_Matrix& A, const BinaryMatrix& C
) const {
    if (C.rows != stat_sec_param || C.cols != ell) {
        std::cout << "plaintext matrix dimension mismatch!\n";
        assert(false);
    }

    std::vector<Multi_Limb_Matrix> sks; 
    sks.reserve(C.rows);
    for (uint64_t i = 0; i < C.rows; i++)
        sks.push_back(preproc_lhe.sampleSecretKey());

    std::vector<Multi_Limb_Matrix> result_cts; result_cts.reserve(C.rows);
    // encrypt each row of C
    for (uint64_t row_ind = 0; row_ind < C.rows; row_ind++) {
        Matrix pt(C.cols, 1);
        for (uint64_t i = 0; i < C.cols; i++)
            pt.data[i] = C.data[row_ind*C.cols + i];

        Multi_Limb_Matrix ct = preproc_lhe.encrypt(A, sks[row_ind], pt);
        result_cts.push_back(ct);
    }

    assert(result_cts.size() == sks.size());

    return make_pair(result_cts, sks);
}

std::pair<std::vector<Matrix>, std::vector<Matrix>>
VeriSimplePIR::PreprocClientMessageNew(
    const Matrix& A, const BinaryMatrix& C
) const {
    if (C.rows != stat_sec_param || C.cols != ell) {
        std::cout << "plaintext matrix dimension mismatch!\n";
        assert(false);
    }

    std::vector<Matrix> sks; 
    sks.reserve(C.rows);
    for (uint64_t i = 0; i < C.rows; i++)
        sks.push_back(lhe.sampleSecretKeyTernary());

    //for (uint64_t i = 0; i < C.rows; i++)
        //sks.push_back(lhe.sampleSecretKey());

    std::vector<Matrix> result_cts; result_cts.reserve(C.rows);
    // encrypt each row of C
    for (uint64_t row_ind = 0; row_ind < C.rows; row_ind++) {
        Matrix pt(C.cols, 1);
        for (uint64_t i = 0; i < C.cols; i++)
            pt.data[i] = C.data[row_ind*C.cols + i];

        Matrix ct = lhe.encrypt(A, sks[row_ind], pt, true);
        result_cts.push_back(ct);
    }

    assert(result_cts.size() == sks.size());

    return make_pair(result_cts, sks);
}

std::pair<std::vector<Matrix64>, std::vector<Matrix64>>
VeriSimplePIR::PreprocClientMessageNew64(
    const Matrix64& A, const BinaryMatrix& C
) const {
    if (C.rows != stat_sec_param || C.cols != ell) {
        std::cout << "plaintext matrix dimension mismatch!\n";
        assert(false);
    }

    std::vector<Matrix64> sks; 
    sks.reserve(C.rows);
    for (uint64_t i = 0; i < C.rows; i++)
        sks.push_back(lhe.sampleSecretKeyTernary64());

    //for (uint64_t i = 0; i < C.rows; i++)
        //sks.push_back(lhe.sampleSecretKey());

    std::vector<Matrix64> result_cts; result_cts.reserve(C.rows);
    // encrypt each row of C
    for (uint64_t row_ind = 0; row_ind < C.rows; row_ind++) {
        Matrix64 pt(C.cols, 1);
        for (uint64_t i = 0; i < C.cols; i++)
            pt.data[i] = C.data[row_ind*C.cols + i];

        Matrix64 ct = lhe.encrypt64(A, sks[row_ind], pt, true);
        result_cts.push_back(ct);
    }

    assert(result_cts.size() == sks.size());

    return make_pair(result_cts, sks);
}

std::pair<std::vector<Multi_Limb_Matrix>, std::vector<Multi_Limb_Matrix>>
VeriSimplePIR::PreprocFakeClientMessage() const {
    std::vector<Multi_Limb_Matrix> sks; 
    sks.reserve(stat_sec_param);
    for (uint64_t i = 0; i < stat_sec_param; i++)
        sks.push_back(preproc_lhe.sampleSecretKey());

    std::vector<Multi_Limb_Matrix> result_cts; result_cts.reserve(sks.size());
    // encrypt each row of C
    for (uint64_t row_ind = 0; row_ind < sks.size(); row_ind++) {
        Multi_Limb_Matrix ct(ell, 1);
        random_fast(ct.q_data); 
        random_fast(ct.kappa_data, preproc_lhe.kappa); 
        result_cts.push_back(ct);
    }

    assert(result_cts.size() == sks.size());

    return make_pair(result_cts, sks);
}

// Just multiply input ciphertexts by D^T
std::vector<Multi_Limb_Matrix> VeriSimplePIR::PreprocAnswer(
    const std::vector<Multi_Limb_Matrix>& in_cts, 
    const Matrix& D
    // const PackedMatrix& D
) const {
    // if (D.orig_rows != m || D.orig_cols != ell) {
    //     std::cout << "plaintext matrix dimension mismatch!\n";
    //     assert(false);
    // }

    std::vector<Multi_Limb_Matrix> result_cts; 
    result_cts.reserve(in_cts.size());
    for (uint64_t i = 0; i < in_cts.size(); i++)
        result_cts.push_back(matMulVec(D, in_cts[i], preproc_lhe.kappa));
        // result_cts.push_back(matVecMulColPacked(D, in_cts[i], preproc_lhe.kappa));

    return result_cts;
}

std::vector<Matrix> VeriSimplePIR::PreprocAnswerNew(
    const std::vector<Matrix>& in_cts, 
    const Matrix& D
    // const PackedMatrix& D
) const {
    // if (D.orig_rows != m || D.orig_cols != ell) {
    //     std::cout << "plaintext matrix dimension mismatch!\n";
    //     assert(false);
    // }

    std::vector<Matrix> result_cts; 
    result_cts.reserve(in_cts.size());
    for (uint64_t i = 0; i < in_cts.size(); i++)
        result_cts.push_back(matMulVec(D, in_cts[i]));
        // result_cts.push_back(matVecMulColPacked(D, in_cts[i], preproc_lhe.kappa));

    return result_cts;
}

// takes in one row of D^T and multiplies it m times
// used when database is too big for benchmarking machine
std::vector<Multi_Limb_Matrix> VeriSimplePIR::PreprocFakeComputeAnswer(
    const std::vector<Multi_Limb_Matrix>& in_cts, 
    const Matrix& D
) const {
    std::vector<Multi_Limb_Matrix> result_cts; 
    result_cts.reserve(in_cts.size());
    for (uint64_t i = 0; i < in_cts.size(); i++) {
        Multi_Limb_Matrix ct_elem(1, 1);
        for (uint64_t j = 0; j < m; j++) {
            ct_elem = matMulVec(D, in_cts[i], preproc_lhe.kappa);
        }
        Multi_Limb_Matrix ct(ell, 1);
        for (uint64_t j = 0; j < ell; j++) {
            ct.q_data.data[j] = ct_elem.q_data.data[0];
            ct.kappa_data.data[j] = ct_elem.kappa_data.data[0];
        }
        result_cts.push_back(ct);
    }

    assert(result_cts.size() == in_cts.size());
    return result_cts;
}


std::vector<Multi_Limb_Matrix> VeriSimplePIR::PreprocFakeAnswer() const {
    std::vector<Multi_Limb_Matrix> result_cts; 
    result_cts.reserve(stat_sec_param);
    for (uint64_t i = 0; i < stat_sec_param; i++) {
        Multi_Limb_Matrix ct(m, 1);
        random_fast(ct.q_data); 
        random_fast(ct.kappa_data, preproc_lhe.kappa); 
        result_cts.push_back(ct);
    }
    return result_cts;
}

void VeriSimplePIR::HashAandH(unsigned char * hash, const Multi_Limb_Matrix& A, const Multi_Limb_Matrix& H) const {
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256, A.q_data.data, A.rows*A.cols*sizeof(Elem));
    SHA256_Update(&sha256, A.kappa_data.data, A.rows*A.cols*sizeof(Elem));
    SHA256_Update(&sha256, H.q_data.data, H.rows*H.cols*sizeof(Elem));
    SHA256_Update(&sha256, H.kappa_data.data, H.rows*H.cols*sizeof(Elem));
    SHA256_Final(hash, &sha256);
}

// This is the C used to prove the correctness of the preprocessed computation. 
// The dimension is lambda x m
BinaryMatrix VeriSimplePIR::BatchHashToC(const unsigned char * AandHhash, const std::vector<Multi_Limb_Matrix>& u_vec, const std::vector<Multi_Limb_Matrix>& v_vec) const {
    
    assert(u_vec.size() == v_vec.size());

    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256_CTX sha256;
    SHA256_Init(&sha256);

    for (uint64_t i = 0; i < u_vec.size(); i++) {
        SHA256_Update(&sha256, u_vec[i].q_data.data, u_vec[i].rows*u_vec[i].cols*sizeof(Elem));
        SHA256_Update(&sha256, u_vec[i].kappa_data.data, u_vec[i].rows*u_vec[i].cols*sizeof(Elem));
        SHA256_Update(&sha256, v_vec[i].q_data.data, v_vec[i].rows*v_vec[i].cols*sizeof(Elem));
        SHA256_Update(&sha256, v_vec[i].kappa_data.data, v_vec[i].rows*v_vec[i].cols*sizeof(Elem));
    }
    SHA256_Final(hash, &sha256);

    SeedType seed = osuCrypto::toBlock(hash) ^ osuCrypto::toBlock(hash + 16);
    // memcpy((unsigned char *)seed, hash, sizeof(seed)); 
    // std::cout << "Dummy C seed!\n";
    BinaryMatrix C(stat_sec_param, m);
    pseudorandom(C, seed);
    return C;
}

// Prove that the computation in PreprocAnswer was correct.
Matrix VeriSimplePIR::PreprocProve(
    const unsigned char * hash,
    const std::vector<Multi_Limb_Matrix>& u, const std::vector<Multi_Limb_Matrix>& v, 
    // const PackedMatrix& D
    const Matrix& D
) const {
    BinaryMatrix C = BatchHashToC(hash, u, v);
    // Matrix Z = matMulLeftBinaryRightColPacked_Hardcoded(C, D);
    Matrix Z = matMulLeftBinary(C, D);
    return Z;
}

Matrix VeriSimplePIR::PreprocFakeProve() const {
    Matrix Z(stat_sec_param, ell);
    random(Z, lhe.p*m);
    return Z;
}

// Verifies preprocessed Z
void VeriSimplePIR::PreprocVerify(
    const Multi_Limb_Matrix& A, const Multi_Limb_Matrix& H, 
    const unsigned char * hash,
    const std::vector<Multi_Limb_Matrix>& u, const std::vector<Multi_Limb_Matrix>& v, 
    const Matrix& Z, const bool fake
) const {
    const size_t norm_bound = lhe.p*m;  // D^T rows have length m
    const size_t Z_len = Z.rows * Z.cols;
    for (size_t i = 0; i < Z_len; i++) {
        if (Z.data[i] >= norm_bound) {
            if (!fake){
                std::cout << "Z is too big!\n";
                assert(false);
            }
        }
    }

    BinaryMatrix C = BatchHashToC(hash, u, v);
    
    Multi_Limb_Matrix leftMatFixed = matMul(Z, A, preproc_lhe.kappa);
    Multi_Limb_Matrix rightMatFixed = matMulLeftBinary(C, H, preproc_lhe.kappa);
    if (!eq(leftMatFixed, rightMatFixed, !fake)) {
        if (!fake){
            std::cout << "preproc verify mismatch!\n";
            assert(false);
        }
    }

    for (uint64_t i = 0; i < u.size(); i++) {
        Multi_Limb_Matrix leftMatQ = matMul(Z, u[i], preproc_lhe.kappa);
        Multi_Limb_Matrix rightMatQ = matMulLeftBinary(C, v[i], preproc_lhe.kappa);
        if (!eq(leftMatQ, rightMatQ, !fake)) {
            if (!fake){
                std::cout << "preproc verify mismatch!\n";
                assert(false);
            }
        }
    }
}

// Decrypts result Z
// Check Z against A_1 and plaintext C
Matrix VeriSimplePIR::PreprocRecoverZ(
    const Multi_Limb_Matrix& H_2,
    const std::vector<Multi_Limb_Matrix>& sks,
    const std::vector<Multi_Limb_Matrix>& res_ct
) const {
    assert(res_ct.size() == stat_sec_param);
    assert(res_ct.size() == sks.size());
    Matrix Z(stat_sec_param, m);

    for (uint64_t row_ind = 0; row_ind < res_ct.size(); row_ind++) {
        Matrix Z_row = preproc_lhe.decrypt(H_2, sks[row_ind], res_ct[row_ind]);
        assert(Z_row.rows == m); assert(Z_row.cols == 1);
        
        for (uint64_t col_ind = 0; col_ind < m; col_ind++)
            Z.data[row_ind*m + col_ind] = Z_row.data[col_ind];
    }

    return Z;
}

Matrix VeriSimplePIR::PreprocRecoverZNew(
    const Matrix& H_2,
    const std::vector<Matrix>& sks,
    const std::vector<Matrix>& res_ct,
    uint32_t blowup
) const {
    assert(res_ct.size() == stat_sec_param);
    assert(res_ct.size() == sks.size());
    //std::cout << "Ct size: " << res_ct.size() << std::endl;
    //std::cout << "Sk size: " << sks.size() << std::endl;
    Matrix Z(stat_sec_param, m);

    for (uint64_t row_ind = 0; row_ind < res_ct.size(); row_ind++) {
        //Matrix Z_row = lhe.decrypt(H_2, sks[row_ind], res_ct[row_ind], true);
        Matrix Z_row = lhe.decrypt_preproc(H_2, sks[row_ind], res_ct[row_ind], blowup, true);
        assert(Z_row.rows == m); assert(Z_row.cols == 1);
        
        for (uint64_t col_ind = 0; col_ind < m; col_ind++)
            Z.data[row_ind*m + col_ind] = Z_row.data[col_ind];
    }

    return Z;
}

Matrix VeriSimplePIR::PreprocRecoverZNew64(
    const Matrix64& H_2,
    const std::vector<Matrix64>& sks,
    const std::vector<Matrix64>& res_ct,
    uint32_t blowup
) const {
    assert(res_ct.size() == stat_sec_param);
    assert(res_ct.size() == sks.size());
    //std::cout << "Ct size: " << res_ct.size() << std::endl;
    //std::cout << "Sk size: " << sks.size() << std::endl;
    Matrix Z(stat_sec_param, m);

    for (uint64_t row_ind = 0; row_ind < res_ct.size(); row_ind++) {
        //Matrix Z_row = lhe.decrypt(H_2, sks[row_ind], res_ct[row_ind], true);
        Matrix Z_row = lhe.decrypt_preproc64(H_2, sks[row_ind], res_ct[row_ind], blowup, true);
        assert(Z_row.rows == m); assert(Z_row.cols == 1);
        
        for (uint64_t col_ind = 0; col_ind < m; col_ind++)
            Z.data[row_ind*m + col_ind] = Z_row.data[col_ind];
    }

    return Z;
}

void VeriSimplePIR::VerifyPreprocZ(
    const Matrix& Z,
    const Matrix& A_1, const BinaryMatrix& C, const Matrix& H_1,
    const bool fake
) const {
    // Verify Z against A_1 and H_1
    // std::cout << "Z = "; print(Z); 

    const size_t norm_bound = lhe.p*ell; 
    const size_t Z_len = Z.rows * Z.cols;
    for (size_t i = 0; i < Z_len; i++) {
        if (Z.data[i] >= norm_bound && !fake) {
            std::cout << "Z is too big!\n";
            assert(false);
        }
    }

    //std::cout << "A = "; print(A_1);

    Matrix leftMat = matMul(Z, A_1);
    Matrix rightMat = matMulLeftBinary(C, H_1);
    if (!eq(leftMat, rightMat) && !fake) {
        std::cout << "\n\nERROR: \n Primary plaintext verify mismatch!\n\n";
        assert(false);
    }
    else {
        std::cout << "Primary verify success\n";
    }
}

void VeriSimplePIR::VerifyPreprocZCompressed(
    const Matrix& RZ,
    const Matrix& A_1, const Matrix& RC, const Matrix& H_1,
    const bool fake
) const {

    //std::cout << "A = "; print(A_1);

    Matrix leftMat = matMul(RZ, A_1);
    Matrix rightMat = matMul(RC, H_1);
    if (!eq(leftMat, rightMat) && !fake) {
        std::cout << "\n\nERROR: \n Primary plaintext verify mismatch!\n\n";
        assert(false);
    }
    else {
        std::cout << "Primary verify success\n";
    }
}

void VeriSimplePIR::VerifyPreprocZSecondary(
    const Matrix& Z,
    const Matrix& A_1, const BinaryMatrix& C, const Matrix& H_1,
    const bool fake
) const { 

    const size_t norm_bound = lhe.p*ell; 
    const size_t Z_len = Z.rows * Z.cols;
    for (size_t i = 0; i < Z_len; i++) {
        if (Z.data[i] >= norm_bound && !fake) {
            std::cout << "Z is too big!\n";
            assert(false);
        }
    }

    Matrix leftMat = matMulSecondary(Z, A_1);
    Matrix rightMat = matMulLeftBinarySecondary(C, H_1);
    if (!eq(leftMat, rightMat) && !fake) {
        std::cout << "\n\nERROR: \n Secondary plaintext verify mismatch!\n\n";
        assert(false);
    }
    else {
        std::cout << "Secondary verify success\n";
    }
}

void VeriSimplePIR::VerifyPreprocZSecondaryCompressed(
    const Matrix& RZ,
    const Matrix& A_1, const Matrix& RC, const Matrix& H_1,
    const bool fake
) const { 

    Matrix leftMat = matMulSecondary(RZ, A_1);
    Matrix rightMat = matMulSecondary(RC, H_1);
    if (!eq(leftMat, rightMat) && !fake) {
        std::cout << "\n\nERROR: \n Secondary plaintext verify mismatch!\n\n";
        assert(false);
    }
    else {
        std::cout << "Secondary verify success\n";
    }
}

// Operations for online phase

Matrix VeriSimplePIR::Init() const {
    return lhe.genPublicA(m);
}

Matrix VeriSimplePIR::EmptyA() const {
    Matrix A(m, lhe.n);
    return A;
}

Matrix VeriSimplePIR::FakeInit() const {
    // return lhe.genPublicA(m);
    Matrix A(m, lhe.n);
    random_fast(A);
    return A;
}

Matrix VeriSimplePIR::GenerateHint(const Matrix& A, const Matrix& D, bool solinas, int record_bitsize) const {
    if (D.rows != ell || D.cols != m) {
        std::cout << "database dimension mismatch!\n";
        assert(false);
    }

    /*
    std::cout << "HERE NOWwwww Rows: " << A.rows << " Cols: " << A.cols << std::endl;
    for (uint64_t i = 0; i < A.rows; i++){
        for (uint64_t j = 0; j < A.cols; j++){
            if (j == 0) std::cout << " " << A.data[i*A.cols + j];
        }
    }
    std::cout << std::endl;
    */

    Matrix H = matMul(D, A);
    return H;
}

Matrix VeriSimplePIR::GenerateFakeHint() const {
    Matrix H(ell, lhe.n);
    random_fast(H);
    return H;
}

std::pair<Matrix, Matrix> VeriSimplePIR::Query(const Matrix& A, const uint64_t index) const {
    if (index >= N) {
        std::cout << "index out of range!\n";
        assert(false);
    }

    const uint64_t index_col = dbParams.indexToColumn(index);
    // std::cout << "Query index column = " << index_col << std::endl;

    Matrix pt(m, 1);
    constant(pt, 0);
    pt.data[index_col] = 1;

    Matrix secretKey = lhe.sampleSecretKey();

    Matrix ciphertext = lhe.encrypt(A, secretKey, pt);

    return std::make_pair(ciphertext, secretKey);
}

Matrix VeriSimplePIR::QueryGivenLweSk(const Matrix& A, const uint64_t index, const Matrix& secretKey_in) const {
    if (index >= N) {
        std::cout << "index out of range!\n";
        assert(false);
    }

    const uint64_t index_col = dbParams.indexToColumn(index);
    std::cout << "Query index column = " << index_col << std::endl;

    Matrix pt(m, 1);
    constant(pt, 0);
    pt.data[index_col] = 1;

    Matrix ciphertext = lhe.encrypt(A, secretKey_in, pt);

    return ciphertext;
}

Matrix VeriSimplePIR::QueryGivenLweSkPNNS(const Matrix& A, const uint64_t index, const Matrix& secretKey_in) const {
    if (index >= N) {
        std::cout << "index out of range!\n";
        assert(false);
    }

    // This function is just a testing function leading up to PNNS development.
    const uint64_t index_col = dbParams.indexToColumnPNNS(index);
    std::cout << "Query index column = " << index_col << std::endl;

    Matrix pt(m, 1);
    constant(pt, 0);
    pt.data[index_col] = 1;

    Matrix ciphertext = lhe.encrypt(A, secretKey_in, pt);

    return ciphertext;
}

std::pair<Matrix, Matrix> VeriSimplePIR::RandomPNNSQueryGivenLweSk(const Matrix& A, const Matrix& secretKey_in) const {

    Matrix pt(m, 1);
    constant(pt, 0);
    Elem mask = (1ULL << d) - 1;
    for (uint64_t i = 0; i < m; i++) {
        pt.data[i] = std::rand() & mask;
    }

    Matrix ciphertext = lhe.encrypt(A, secretKey_in, pt);

    return std::pair(ciphertext, pt);
}

Matrix VeriSimplePIR::Answer(const Matrix& ciphertext, const Matrix& D) const {
    if (ciphertext.cols == 1) {
        Matrix ans = matMulVec(D, ciphertext);
        return ans;
    } else {
        Matrix ans = matMul(D, ciphertext);
        return ans;
    }
}

Matrix VeriSimplePIR::Answer(const Matrix& ciphertext, const PackedMatrix& D_packed) const {
    if (ciphertext.cols == 1) {
        Matrix ans = simplepir_matVecMulColPacked_variableCompression(D_packed, ciphertext);
        return ans;
    } else {
        Matrix ans = matMulColPacked(D_packed, ciphertext);
        return ans;
    }
}

void VeriSimplePIR::PreVerify(const Matrix& u, const Matrix& v, const Matrix& Z, const BinaryMatrix& C, const bool fake) const {
    const auto left = matMulVec(Z, u);
    //std::cout << "Done left\n";
    const auto right = matBinaryMulVec(C, v);
    //std::cout << "Done right\n";
    if (!eq(left, right, true)) {
        if (!fake) {
            std::cout << "verify mismatch!\n";
            assert(false);
        }
    }
    //std::cout << "Done PreVerify" << std::endl;
    return;
}

void VeriSimplePIR::PreVerifyCompressed(const Matrix& u, const Matrix& v, const Matrix& rZ, const Matrix& rC, const bool fake) const {
    const auto left = matMulVec(rZ, u);
    //std::cout << "Done left\n";
    const auto right = matMulVec(rC, v);
    //std::cout << "Done right\n";
    if (!eq(left, right, true)) {
        std::cout << "\n\nERROR:\n D*u verify mismatch!\n";
        if (!fake) {
            std::cout << "\n\nERROR:\n D*u verify mismatch!\n";
            assert(false);
        }
    }
    else {
        std::cout << "D*u verify success\n";
    }
    //std::cout << "Done Compressed PreVerify" << std::endl;
    return;
}

void VeriSimplePIR::HsVerify(const Matrix& s, const Matrix& w, const Matrix& Z, const BinaryMatrix& C, const Matrix& A, const bool fake) const {
    const auto As = matMulVec(A, s);
    const auto left = matMulVec(Z, As);
    //std::cout << "Done left\n";
    const auto right = matBinaryMulVec(C, w);
    //std::cout << "Done right\n";
    if (!eq(left, right, true)) {
        if (!fake) {
            std::cout << "verify mismatch!\n";
            assert(false);
        }
    }
    //std::cout << "Done our new H*s verification" << std::endl;
    return;
}

void VeriSimplePIR::HsVerifyCompressed(const Matrix& s, const Matrix& w, const Matrix& rZ, const Matrix& rC, const Matrix& A, const bool fake) const {
    const auto As = matMulVec(A, s);
    const auto left = matMulVec(rZ, As);
    //std::cout << "Done left\n";
    const auto right = matMulVec(rC, w);
    //std::cout << "Done right\n";
    if (!eq(left, right, true)) {
        std::cout << "\n\nERROR:\n H*s verify mismatch!\n";
        if (!fake) {
            std::cout << "\n\nERROR:\n H*s verify mismatch!\n";
            assert(false);
        }
    }
    else {
            std::cout << "H*s verify success\n";
    }
    //std::cout << "Done our new compressed H*s verification" << std::endl;
    return;
}

void VeriSimplePIR::HsVerifyCompressed(const Matrix& As, const Matrix& w, const Matrix& rZ, const Matrix& rC, const bool fake) const {
    const auto left = matMulVec(rZ, As);
    //std::cout << "Done left\n";
    const auto right = matMulVec(rC, w);
    //std::cout << "Done right\n";
    if (!eq(left, right, true)) {
        std::cout << "\n\nERROR:\n H*s verify mismatch!\n";
        if (!fake) {
            std::cout << "\n\nERROR:\n H*s verify mismatch!\n";
            assert(false);
        }
    }
    else {
            std::cout << "H*s verify success\n";
    }
    //std::cout << "Done our new compressed H*s verification" << std::endl;
    return;
}

void VeriSimplePIR::FakePreVerify(const Matrix& u, const Matrix& v, const Matrix& Z, const BinaryMatrix& C) const {
    PreVerify(u, v, Z, C, true);
}

std::pair<BinaryMatrix, Matrix> VeriSimplePIR::SampleFakeCandZ() const {
    // C is binary with dimension lambda x ell
    // Z is random lambda x m with inf norm at most p*ell
    BinaryMatrix C(stat_sec_param, ell); random(C);
    Matrix Z(stat_sec_param, m); random(Z, lhe.p*ell);

    return std::make_pair(C, Z);
}

entry_t VeriSimplePIR::Recover(const Matrix& hint, const Matrix& ciphertext, const Matrix& secretKey, const uint64_t index) const {
    // const uint64_t index_row = index / m;
    const uint64_t index_row = dbParams.indexToRow(index);
    std::cout << "Recover index row = " << index_row << std::endl;
    if (index_row >= ell) {
        std::cout << "index row is too big!\n";
        assert(false);
    }
  
    Matrix pt = lhe.decrypt(hint, secretKey, ciphertext);
    std::cout << "Recover pt =\n"; print(pt);
    return dbParams.recover(&pt.data[index_row], index);
}

entry_t VeriSimplePIR::RecoverPNNS(const Matrix& hint, const Matrix& ciphertext, const Matrix& secretKey, const uint64_t index) const {
    // const uint64_t index_row = index / m;
    const uint64_t index_row = dbParams.indexToRowPNNS(index);
    std::cout << "Recover index row = " << index_row << std::endl;
    if (index_row >= ell) {
        std::cout << "index row is too big!\n";
        assert(false);
    }
  
    Matrix pt = lhe.decrypt(hint, secretKey, ciphertext);
    //std::cout << "Recover pt =\n"; print(pt);
    return dbParams.recover(&pt.data[index_row], index);
}

Matrix VeriSimplePIR::RecoverPNNSScores(const Matrix& hint, const Matrix& ciphertext, const Matrix& secretKey) const {
    Matrix pt = lhe.decrypt(hint, secretKey, ciphertext);
    //std::cout << "Recover pt =\n"; print(pt);
    return pt;
}