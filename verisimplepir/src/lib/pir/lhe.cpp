#include "lhe.h"

// input is number of rows
Matrix LHE::genPublicA(uint64_t m) const {
    Matrix A(m, n);
    random(A);
    for (uint64_t i = 0; i < A.rows; i++){
        for (uint64_t j = 0; j < A.cols; j++){
            A.data[i*A.cols + j] = fast_mod_max_64(A.data[i*A.cols + j]);
        }
    }
    return A;
};

// input is number of rows
Matrix LHE::genPublicA(uint64_t m, uint64_t n_) const {
    Matrix A(m, n_);
    random(A);
    for (uint64_t i = 0; i < A.rows; i++){
        for (uint64_t j = 0; j < A.cols; j++){
            A.data[i*A.cols + j] = fast_mod_max_64(A.data[i*A.cols + j]);
        }
    }
    return A;
};

// input is number of rows
Matrix64 LHE::genPublicA64(uint64_t m) const {
    Matrix64 A(m, n);
    random64(A);
    for (uint64_t i = 0; i < A.rows; i++){
        for (uint64_t j = 0; j < A.cols; j++){
            A.data[i*A.cols + j] = fast_mod_vspir_max_64(A.data[i*A.cols + j]);
        }
    }
    return A;
};

// input is number of rows
Matrix LHE::genPublicASecondary(uint64_t m) const {
    Matrix A(m, n);
    random(A);
    for (uint64_t i = 0; i < A.rows; i++){
        for (uint64_t j = 0; j < A.cols; j++){
            A.data[i*A.cols + j] = fast_mod_secondary_max_64(A.data[i*A.cols + j]);
        }
    }
    return A;
};

// input is number of rows
Matrix LHE::genPublicASecondary(uint64_t m, uint64_t n_) const {
    Matrix A(m, n_);
    random(A);
    for (uint64_t i = 0; i < A.rows; i++){
        for (uint64_t j = 0; j < A.cols; j++){
            A.data[i*A.cols + j] = fast_mod_secondary_max_64(A.data[i*A.cols + j]);
        }
    }
    return A;
};

// column vector where # of rows is always n
// TODO: need rejection sampling here to ensure uniformity. Not called in our current code
Matrix LHE::sampleSecretKey() const {
    Matrix sk(n, 1);
    random(sk);
    for (uint64_t i = 0; i < sk.rows; i++){
        for (uint64_t j = 0; j < sk.cols; j++){
            //std::cout << sk.data[i*sk.cols + j] << " ";
            sk.data[i*sk.cols + j] = fast_mod_max_64(sk.data[i*sk.cols + j]);
            //std::cout << sk.data[i*sk.cols + j] << std::endl;

        }
    }
    return sk;
};

// column vector where # of rows is always n
Matrix LHE::sampleSecretKeyTernary() const {
    Matrix sk(n, 1);
    randomTernary(sk);
    Elem two = 2;
    for (uint64_t i = 0; i < sk.rows; i++){
        for (uint64_t j = 0; j < sk.cols; j++){
            //std::cout << sk.data[i*sk.cols + j] << " ";
            //sk.data[i*sk.cols + j] = fast_mod_max_64(sk.data[i*sk.cols + j]);
            if (sk.data[i*sk.cols + j] == two) {
                sk.data[i*sk.cols + j] = (Elem)SOL_NTT_PRIME - 1;
            }
            //std::cout << sk.data[i*sk.cols + j] << std::endl;
        }
    }
    return sk;
};

// column vector where # of rows is always n
Matrix64 LHE::sampleSecretKeyTernary64() const {
    Matrix64 sk(n, 1);
    randomTernary64(sk);
    Elem64 two = 2;
    for (uint64_t i = 0; i < sk.rows; i++){
        for (uint64_t j = 0; j < sk.cols; j++){
            //std::cout << sk.data[i*sk.cols + j] << " ";
            //sk.data[i*sk.cols + j] = fast_mod_max_64(sk.data[i*sk.cols + j]);
            if (sk.data[i*sk.cols + j] == two) {
                sk.data[i*sk.cols + j] = (Elem64)VSPIR_LHE_PRIME - 1;
            }
            //std::cout << sk.data[i*sk.cols + j] << std::endl;
        }
    }
    return sk;
};

void LHE::randomPlaintext(Matrix& pt) const {
    random(pt, p);
}

// plaintext has length m, where m is in the parameter used to sample m
// Assuming all elements of pt are less than p;
Matrix LHE::encrypt(const Matrix& A, const Matrix& sk, const Matrix& pt, bool preproc) const {
    if (A.rows != pt.rows) {
        std::cout << "Plaintext dimension mismatch!\n";
        assert(false);
    }

    if (sk.cols != 1 || pt.cols != 1) {
        std::cout << "secret key or plaintext are not column vectors!\n";
        assert(false);
    }

    Matrix ciphertext(A.rows, 1);
    error(ciphertext);
    // constant(&ciphertext); std::cout << "change me back!\n";


    Matrix A_sk = matMulVec(A, sk); 

    matAddInPlace(ciphertext, A_sk);

    // scale plaintext
    Elem Delta__ = 0;
    if (preproc) {
        Delta__ = Delta_preproc;
    }
    else{
        Delta__ = Delta;
    }
    Matrix pt_scaled = matMulScalar(pt, Delta__);

    matAddInPlace(ciphertext, pt_scaled);

    return ciphertext;
};

// plaintext has length m, where m is in the parameter used to sample m
// Assuming all elements of pt are less than p;
Matrix64 LHE::encrypt64(const Matrix64& A, const Matrix64& sk, const Matrix64& pt, bool preproc) const {
    if (A.rows != pt.rows) {
        std::cout << "Plaintext dimension mismatch!\n";
        assert(false);
    }

    if (sk.cols != 1 || pt.cols != 1) {
        std::cout << "secret key or plaintext are not column vectors!\n";
        assert(false);
    }

    Matrix64 ciphertext(A.rows, 1);
    error64(ciphertext);
    // constant(&ciphertext); std::cout << "change me back!\n";


    Matrix64 A_sk = matMulVec64(A, sk); 

    matAddInPlace64(ciphertext, A_sk);

    // scale plaintext
    Elem64 Delta64__ = 0;
    if (preproc) {
        Delta64__ = Delta64_preproc;
    }
    else{
        Delta64__ = Delta64;
    }
    Matrix64 pt_scaled = matMulScalar64(pt, Delta64__);

    matAddInPlace64(ciphertext, pt_scaled);

    return ciphertext;
};

// length of ct should match the # of rows of H
Matrix LHE::decrypt(const Matrix& H, const Matrix& sk, const Matrix& ct, bool preproc) const {
    if (H.rows != ct.rows) {
        std::cout << "Ciphertext dimension mismatch!\n";
        std::cout << H.rows << " x " << H.cols << " vs. " << ct.rows << " x " << ct.cols << std::endl;
        assert(false);
    }

    if (sk.cols != 1 || ct.cols != 1) {
        std::cout << "secret key or ciphertext are not column vectors!\n";
        assert(false);
    }

    Matrix H_sk = matMulVec(H, sk);

    Matrix scaled_pt = matSub(ct, H_sk);

    Elem Delta__ = 0;
    if (preproc) {
        Delta__ = Delta_preproc;
    }
    else{
        Delta__ = Delta;
    }
    Matrix pt = matDivScalar(scaled_pt, Delta__);
    //Matrix pt = matDivScalar(ct, Delta__); // remove this when changing

    for (size_t i = 0; i < pt.rows; i++) {
        if (pt.data[i] == p) {
            pt.data[i] = 0;
        }
    }
    


    return pt;
}

// length of ct should match the # of rows of H
Matrix LHE::decrypt_preproc(const Matrix& H, const Matrix& sk, const Matrix& ct, uint32_t blowup, bool preproc) const {
    if (H.rows != ct.rows) {
        std::cout << "Ciphertext dimension mismatch!\n";
        std::cout << H.rows << " x " << H.cols << " vs. " << ct.rows << " x " << ct.cols << std::endl;
        assert(false);
    }

    if (sk.cols != 1 || ct.cols != 1) {
        std::cout << "secret key or ciphertext are not column vectors!\n";
        assert(false);
    }

    Matrix H_sk = matMulVec(H, sk);

    Matrix scaled_pt = matSub(ct, H_sk);

    Elem Delta__ = 0;
    if (preproc) {
        Delta__ = Delta_preproc;
    }
    else{
        Delta__ = Delta;
    }
    //Matrix pt = matMultThenDivScalar(scaled_pt, q, p*blowup);
    Matrix pt = matDivScalar(scaled_pt, Delta__);

    for (size_t i = 0; i < pt.rows; i++) {
        if (pt.data[i] == p*blowup) {
            pt.data[i] = 0;
        }
    }

    return pt;
}

// length of ct should match the # of rows of H
Matrix LHE::decrypt_preproc64(const Matrix64& H, const Matrix64& sk, const Matrix64& ct, uint32_t blowup, bool preproc) const {
    if (H.rows != ct.rows) {
        std::cout << "Ciphertext dimension mismatch!\n";
        std::cout << H.rows << " x " << H.cols << " vs. " << ct.rows << " x " << ct.cols << std::endl;
        assert(false);
    }

    if (sk.cols != 1 || ct.cols != 1) {
        std::cout << "secret key or ciphertext are not column vectors!\n";
        assert(false);
    }

    Matrix64 H_sk = matMulVec64(H, sk);

    Matrix64 scaled_pt = matSub64(ct, H_sk);

    Elem64 Delta64__ = 0;
    if (preproc) {
        Delta64__ = Delta64_preproc;
    }
    else{
        Delta64__ = Delta64;
    }
    //Matrix pt = matMultThenDivScalar(scaled_pt, q, p*blowup);
    Matrix pt = matDivScalar64(scaled_pt, Delta64__);

    for (size_t i = 0; i < pt.rows; i++) {
        if (pt.data[i] == p*blowup) {
            pt.data[i] = 0;
        }
    }

    return pt;
}