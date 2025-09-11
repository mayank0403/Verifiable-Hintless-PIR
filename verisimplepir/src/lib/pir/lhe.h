#pragma once
#include "mat.h"

class LHE {  // linearly homomorphic encryption
public:
    // Updating based on new lattice attacks [see HintlessPIR]
    // Change: making this updatable to be able to run smaller experiments for testing
    //static constexpr uint64_t n = (sizeof(Elem) <= 4) ? 1408 : 2048;  // number of columns of A. morally, the security parameter
    static constexpr uint64_t n = 2816; //1408;
    // static constexpr uint64_t n = (sizeof(Elem) <= 4) ? 1024 : 2048;  // number of columns of A. morally, the security parameter
    Elem p;  // plaintext modulus
    // the ciphertext modulus q is implicit in the data type Elem
    //      logq is either 32 or 64 depending on the size of Elem.
    //static constexpr uint64_t logq = 8*sizeof(Elem);
    static constexpr uint64_t logq = LOG_Q;
    static constexpr uint64_t q = SOL_NTT_PRIME;
    static constexpr uint64_t vspir_q = VSPIR_LHE_PRIME;
    Elem Delta;  // q/p
    Elem Delta_preproc; // q/(\ell*p)
    Elem64 Delta64;
    Elem64 Delta64_preproc;

    LHE () {};

    // LHE(uint64_t n_in, Elem p_in) : n(n_in), p(p_in) {
    LHE(const Elem p_in, const Elem r_in = 0) : p(p_in) {
        assert(p_in >= 2);
        if (SOL_NTT_PRIME_BITS == 0) {
            Elem q_half = 1ULL << (logq - 1);
            Delta = q_half/p;
            Delta <<= 1;
            Delta_preproc = q_half/(p*r_in);
            Delta_preproc <<= 1;

            Elem64 vspir_q_half = 1ULL << (VSPIR_LHE_LOG_Q - 1);
            Delta64 = vspir_q_half/p;
            Delta64 <<= 1;
            Delta64_preproc = vspir_q_half/(p*r_in);
            Delta64_preproc <<= 1;
        }
        else {
            Delta = q/p;
            Delta_preproc = q/(p*r_in);

            Delta64 = vspir_q/p;
            Delta64_preproc = vspir_q/(p*r_in);
        }
    };

    Matrix genPublicA(uint64_t m) const;  // input is number of rows
    Matrix genPublicA(uint64_t m, uint64_t n_) const;
    Matrix64 genPublicA64(uint64_t m) const;  // input is number of rows
    Matrix genPublicASecondary(uint64_t m) const;  // input is number of rows
    Matrix genPublicASecondary(uint64_t m, uint64_t n_) const;
    Matrix sampleSecretKey() const;  // column vector where # of rows is always n
    Matrix sampleSecretKeyTernary() const;  // column vector where # of rows is always n
    Matrix64 sampleSecretKeyTernary64() const;  // column vector where # of rows is always n

    void randomPlaintext(Matrix& pt) const;

    // plaintext has length m, where m is in the parameter used to sample m
    Matrix encrypt(const Matrix& A, const Matrix& sk, const Matrix& pt, bool preproc = false) const;  

    // plaintext has length m, where m is in the parameter used to sample m
    Matrix64 encrypt64(const Matrix64& A, const Matrix64& sk, const Matrix64& pt, bool preproc = false) const;  

    // length of ct should match the # of rows of H
    Matrix decrypt(const Matrix& H, const Matrix& sk, const Matrix& ct, bool preproc = false) const;

    Matrix decrypt_preproc(const Matrix& H, const Matrix& sk, const Matrix& ct, uint32_t blowup, bool preproc) const;

    Matrix decrypt_preproc64(const Matrix64& H, const Matrix64& sk, const Matrix64& ct, uint32_t blowup, bool preproc) const;
};
