/*
    Basic matrix definitions and operations
*/
#pragma once

#include <cstdio>
#include <cstdlib>
#include "utils.h"
#include "math/prng.h"

// typedef std::mt19937_64 PRNG;
typedef osuCrypto::PRNG PRNG;
// typedef PRNG::result_type SeedType;
typedef osuCrypto::block SeedType;

#ifdef SH_RUN
typedef uint32_t Elem; // for PIR
typedef uint64_t Elem64; // To ensure that all computation in semi-honest happens over 32-bit always (never anything over 64 bits)
typedef uint64_t ElemDouble; // same reason as above
#else
/* These are the only tested values of Elem */ 
typedef uint32_t Elem; // for PIR
typedef uint64_t Elem64; // THIS STAYS FIXED; USED BY Highway and it is always enough to allow lazy multiplication because one operand is the database itself
// TODO: To support > 32-bit primes, we need keep Elem64 the same, but introduce a new type "Elem_double" which always is double the space of Elem. And all multiplications where both operands are Elem elements, we use Elem_double (except in inner_product_highway)
typedef uint64_t ElemDouble; // for multiplication of two Z_q values
#endif // SH_RUN

#define ALIGN 8


class Matrix {
public:
    uint64_t rows, cols;
    Elem* data;  // packed in row-major order by default

    Matrix() {};

    // used on unitialized matrices to write data directly 
    void init_no_memset(uint64_t r, uint64_t c) {
        rows = r;
        cols = c;

        #ifdef ALIGN
        data = (Elem*)aligned_alloc(ALIGN, rows*cols * sizeof(Elem));
        #else
        data = (Elem*)malloc(rows*cols * sizeof(Elem));
        #endif
    }

    Matrix(uint64_t r, uint64_t c, const Elem val = 0) : rows(r), cols(c) {
         #ifdef ALIGN
        data = (Elem*)aligned_alloc(ALIGN, rows*cols * sizeof(Elem));
        #else
        data = (Elem*)malloc(rows*cols * sizeof(Elem));
        #endif

        memset(data, val, rows*cols * sizeof(Elem));
    }

    // ~Matrix() {
    //     free(data);
    // }

    Matrix(const Matrix& rhs) {
        rows = rhs.rows;
        cols = rhs.cols;

        #ifdef ALIGN
        data = (Elem*)aligned_alloc(ALIGN, rows*cols * sizeof(Elem));
        #else
        data = (Elem*)malloc(rows*cols * sizeof(Elem));
        #endif
        
        memcpy(data, rhs.data, rows*cols * sizeof(Elem));
    }

    Elem getElem(const uint64_t row, const uint64_t col) const {
        // assumes row-major order packing
        if (row >= rows) {
            std::cout << "row too big!\n";
            std::cout << row << " " << rows << std::endl;
            assert(false);
        }
        if (col > cols) {
            std::cout << "col too big!\n";
            std::cout << col << " " << cols << std::endl;
            assert(false);
        }
        return data[row*cols + col];
    }

    Matrix getColumn(const uint64_t col) const {
        Matrix res(rows, 1);
        for (uint64_t i = 0; i < rows; i++) 
            res.data[i] = data[i*cols + col];
        return res;
    }
};

class Matrix64 {
public:
    uint64_t rows, cols;
    Elem64* data;  // packed in row-major order by default

    Matrix64() {};

    // used on unitialized matrices to write data directly 
    void init_no_memset(uint64_t r, uint64_t c) {
        rows = r;
        cols = c;

        #ifdef ALIGN
        data = (Elem64*)aligned_alloc(ALIGN, rows*cols * sizeof(Elem64));
        #else
        data = (Elem64*)malloc(rows*cols * sizeof(Elem64));
        #endif
    }

    Matrix64(uint64_t r, uint64_t c, const Elem64 val = 0) : rows(r), cols(c) {
         #ifdef ALIGN
        data = (Elem64*)aligned_alloc(ALIGN, rows*cols * sizeof(Elem64));
        #else
        data = (Elem64*)malloc(rows*cols * sizeof(Elem64));
        #endif

        memset(data, val, rows*cols * sizeof(Elem64));
    }

    // ~Matrix64() {
    //     free(data);
    // }

    Matrix64(const Matrix64& rhs) {
        rows = rhs.rows;
        cols = rhs.cols;

        #ifdef ALIGN
        data = (Elem64*)aligned_alloc(ALIGN, rows*cols * sizeof(Elem64));
        #else
        data = (Elem64*)malloc(rows*cols * sizeof(Elem64));
        #endif
        
        memcpy(data, rhs.data, rows*cols * sizeof(Elem64));
    }

    Elem64 getElem(const uint64_t row, const uint64_t col) const {
        // assumes row-major order packing
        if (row >= rows) {
            std::cout << "row too big!\n";
            std::cout << row << " " << rows << std::endl;
            assert(false);
        }
        if (col > cols) {
            std::cout << "col too big!\n";
            std::cout << col << " " << cols << std::endl;
            assert(false);
        }
        return data[row*cols + col];
    }

    Matrix64 getColumn(const uint64_t col) const {
        Matrix64 res(rows, 1);
        for (uint64_t i = 0; i < rows; i++) 
            res.data[i] = data[i*cols + col];
        return res;
    }
};

class MatrixDouble {
public:
    uint64_t rows, cols;
    ElemDouble* data;  // packed in row-major order by default

    MatrixDouble() {};

    // used on unitialized matrices to write data directly 
    void init_no_memset(uint64_t r, uint64_t c) {
        rows = r;
        cols = c;

        #ifdef ALIGN
        data = (ElemDouble*)aligned_alloc(ALIGN, rows*cols * sizeof(ElemDouble));
        #else
        data = (ElemDouble*)malloc(rows*cols * sizeof(ElemDouble));
        #endif
    }

    MatrixDouble(uint64_t r, uint64_t c, const ElemDouble val = 0) : rows(r), cols(c) {
         #ifdef ALIGN
        data = (ElemDouble*)aligned_alloc(ALIGN, rows*cols * sizeof(ElemDouble));
        #else
        data = (ElemDouble*)malloc(rows*cols * sizeof(ElemDouble));
        #endif

        std::cout << "\nNOT SUPPORTED\n" << std::endl;
        assert(false);
        memset(data, 0, rows*cols * sizeof(ElemDouble));
        //memset(data, val, rows*cols * sizeof(ElemDouble));
    }

    MatrixDouble(const MatrixDouble& rhs) {
        rows = rhs.rows;
        cols = rhs.cols;

        #ifdef ALIGN
        data = (ElemDouble*)aligned_alloc(ALIGN, rows*cols * sizeof(ElemDouble));
        #else
        data = (ElemDouble*)malloc(rows*cols * sizeof(ElemDouble));
        #endif
        
        memcpy(data, rhs.data, rows*cols * sizeof(ElemDouble));
    }

    ElemDouble getElem(const uint64_t row, const uint64_t col) const {
        // assumes row-major order packing
        if (row >= rows) {
            std::cout << "row too big!\n";
            std::cout << row << " " << rows << std::endl;
            assert(false);
        }
        if (col > cols) {
            std::cout << "col too big!\n";
            std::cout << col << " " << cols << std::endl;
            assert(false);
        }
        return data[row*cols + col];
    }

    MatrixDouble getColumn(const uint64_t col) const {
        MatrixDouble res(rows, 1);
        for (uint64_t i = 0; i < rows; i++) 
            res.data[i] = data[i*cols + col];
        return res;
    }
};

class BinaryMatrix {
public:
    uint64_t rows, cols;
    bool* data;  // packed in row-major order by default

    BinaryMatrix() {};

    BinaryMatrix(uint64_t r, uint64_t c, const Elem val = 0) : rows(r), cols(c) {
        #ifdef ALIGN
        data = (bool*)aligned_alloc(ALIGN, rows*cols * sizeof(bool));
        #else
        data = (bool*)malloc(rows*cols * sizeof(bool));
        #endif

        memset(data, val, rows*cols * sizeof(bool));
    }

    ~BinaryMatrix() {
        free(data);
    }

    BinaryMatrix(const BinaryMatrix& rhs) {
        rows = rhs.rows;
        cols = rhs.cols;

        #ifdef ALIGN
        data = (bool*)aligned_alloc(ALIGN, rows*cols * sizeof(bool));
        #else
        data = (bool*)malloc(rows*cols * sizeof(bool));
        #endif
        
        memcpy(data, rhs.data, rows*cols * sizeof(bool));
    }

    // used on unitialized matrices to write data directly 
    void init_no_memset(uint64_t r, uint64_t c) {
        rows = r;
        cols = c;

        #ifdef ALIGN
        data = (bool*)aligned_alloc(ALIGN, rows*cols * sizeof(bool));
        #else
        data = (bool*)malloc(rows*cols * sizeof(bool));
        //data = (Elem*)malloc(rows*cols * sizeof(Elem));
        #endif
    }


    Matrix asMatrix() const {
        Matrix res; res.init_no_memset(rows, cols);
        for (size_t i = 0; i < rows*cols; i++)
            res.data[i] = (Elem)data[i];
        return res;
    }

    // Matrix& operator&() {
    //     return (Matrix&)(*this);
    // }
};

void random(Matrix& mat, const Elem max = 0);
void randomTernary(Matrix& mat);
void randomTernary64(Matrix64& mat);
void random64(Matrix64& mat, const Elem64 max = 0);
void random(BinaryMatrix& mat);
void random_fast(Matrix& mat, const Elem modulus = 0);

void pseudorandom(Matrix& mat, const SeedType& seed, const Elem max = 0);
void pseudorandom(BinaryMatrix& mat, const SeedType& seed);

void error(Matrix& mat);
void error64(Matrix64& mat);

template <typename MatrixType>
void constant(MatrixType& mat, const Elem val = 0);

bool eq(const Matrix& a, const Matrix& b, const bool verbose = false);

template <typename MatType>
void print(const MatType& mat);

template <typename MatrixType>
MatrixType transpose(const MatrixType& in) {
    const size_t rows = in.rows;
    const size_t cols = in.cols;

    MatrixType out; out.init_no_memset(cols, rows);

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            out.data[j * rows + i] = in.data[i * cols + j];
        }
    }

    return out;
}

template BinaryMatrix transpose<BinaryMatrix>(const BinaryMatrix& in);
template Matrix transpose<Matrix>(const Matrix& in);


Matrix matAdd(const Matrix &a, const Matrix &b);
void matAddInPlace(Matrix& a, const Matrix& b, const Elem modulus = 0);
void matAddInPlace64(Matrix64& a, const Matrix64& b, const Elem64 modulus = 0);

Matrix matSub(const Matrix &a, const Matrix &b, const Elem modulus = 0);
Matrix64 matSub64(const Matrix64 &a, const Matrix64 &b, const Elem64 modulus = 0);
// void matSubInPlace(Matrix& a, const Matrix& b);

Matrix matMul(const Matrix &a, const Matrix &b, const Elem modulus = 0);
Matrix matMulSecondary(const Matrix &a, const Matrix &b, const Elem modulus = 0);
Matrix matMulLeftBinary(const BinaryMatrix& binary, const Matrix& b, const Elem modulus = 0);
Matrix matMulLeftBinarySecondary(const BinaryMatrix& binary, const Matrix& b, const Elem modulus = 0);
Matrix matMulRightBinary(const Matrix& b, const BinaryMatrix& binary, const Elem modulus = 0);
Matrix matMulRightBinarySecondary(const Matrix& b, const BinaryMatrix& binary, const Elem modulus = 0);


// produces the product a*[b | c], where c is a vector
Matrix matMulAppendVec(const Matrix &a, const Matrix &b, const Matrix &c);
Matrix binaryMatMulAppendVec(const BinaryMatrix &binary, const Matrix &b, const Matrix &c);

Matrix matMulScalar(const Matrix &a, const Elem b, const Elem modulus = 0);
Matrix64 matMulScalar64(const Matrix64 &a, const Elem64 b, const Elem64 modulus = 0);
Matrix matDivScalar(const Matrix &a, const Elem b);
Matrix matDivScalar64(const Matrix64 &a, const Elem64 b);
Matrix matMultThenDivScalar(const Matrix& a, const Elem b, const Elem c);

Matrix matMulVec(const Matrix &a, const Matrix &b, const Elem modulus = 0);
Matrix64 matMulVec64(const Matrix64 &a, const Matrix64 &b, const Elem64 modulus = 0);
Matrix matBinaryMulVec(const BinaryMatrix& a, const Matrix& b);

Matrix RowVecMulMat(const Matrix& v, const Matrix& m);
Matrix RowVecMulBinaryMat(const Matrix& v, const BinaryMatrix& m);