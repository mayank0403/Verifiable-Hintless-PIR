#include "mat.h"
#include "gauss.h"

#include <fstream>
#include <iostream>


void randomInner(Matrix& mat, PRNG& prng, const Elem max) {
    std::uniform_int_distribution<Elem> dist;
    if (max == 0) {
        dist = std::uniform_int_distribution<Elem>(0);  // default max is maximum value of type
    } else {
        dist = std::uniform_int_distribution<Elem>(0, max-1); 
    }
    const size_t len = mat.rows * mat.cols;
    for (size_t i = 0; i < len; i++) {
        mat.data[i] = dist(prng);
    }
}

void randomInner64(Matrix64& mat, PRNG& prng, const Elem64 max) {
    std::uniform_int_distribution<Elem64> dist;
    if (max == 0) {
        dist = std::uniform_int_distribution<Elem64>(0);  // default max is maximum value of type
    } else {
        dist = std::uniform_int_distribution<Elem64>(0, max-1); 
    }
    const size_t len = mat.rows * mat.cols;
    for (size_t i = 0; i < len; i++) {
        mat.data[i] = dist(prng);
    }
}

void randomInner(BinaryMatrix& mat, PRNG& prng) {
    auto dist = std::uniform_int_distribution<>(0, 1); 
    const size_t len = mat.rows * mat.cols;
    for (size_t i = 0; i < len; i++) {
        mat.data[i] = dist(prng);
    }
}

void random_fast(Matrix& mat, const Elem modulus) {
    auto fd = fopen("/dev/random", "r");
    // memcpy(fd, mat.data, mat.rows*mat.cols*sizeof(Elem));
    fgets((char *)mat.data, mat.rows*mat.cols*sizeof(Elem), fd);
    fclose(fd);

    if (modulus > 0) {
        for (size_t i = 0; i < mat.rows*mat.cols; i++)
            mat.data[i] %= modulus;
    }
}

void random(Matrix& mat, const Elem max) {
    PRNG prng(osuCrypto::sysRandomSeed());
    randomInner(mat, prng, max);
}

void randomTernary(Matrix& mat) {
    PRNG prng(osuCrypto::sysRandomSeed());
    randomInner(mat, prng, 3);
}

void randomTernary64(Matrix64& mat) {
    PRNG prng(osuCrypto::sysRandomSeed());
    randomInner64(mat, prng, 3);
}

void random64(Matrix64& mat, const Elem64 max) {
    PRNG prng(osuCrypto::sysRandomSeed());
    randomInner64(mat, prng, max);
}

void random(BinaryMatrix& mat) {
    PRNG prng(osuCrypto::sysRandomSeed());
    randomInner(mat, prng);
}

void pseudorandom(Matrix& mat, const SeedType& seed, const Elem max) {
    PRNG prng(seed);
    randomInner(mat, prng, max);
}

void pseudorandom(BinaryMatrix& mat, const SeedType& seed) {
    PRNG prng(seed);
    randomInner(mat, prng);
}

void error(Matrix& mat) {
    const size_t len = mat.rows * mat.cols;
    //auto error_norm = 0;
    Elem mid_pt = 1 << (sizeof(Elem)*8 - 1);
    Elem magnitude;
    for (size_t i = 0; i < len; i++) {
        mat.data[i] = GaussSample<Elem>();
        if (mat.data[i] >= mid_pt) {
            // Negative value
            magnitude = mat.data[i] * (-1);
            mat.data[i] = (Elem)SOL_NTT_PRIME - magnitude;
        }
        //if ((Elem)mat.data[i] > error_norm) {
        //    error_norm = (Elem)mat.data[i];
        //}
    }
    //std::cout << "error norm = " << error_norm << std::endl;
}

void error64(Matrix64& mat) {
    const size_t len = mat.rows * mat.cols;
    //auto error_norm = 0;
    Elem64 mid_pt = 1 << (sizeof(Elem64)*8 - 1);
    Elem64 magnitude;
    for (size_t i = 0; i < len; i++) {
        mat.data[i] = GaussSample<Elem64>();
        if (mat.data[i] >= mid_pt) {
            // Negative value
            magnitude = mat.data[i] * (-1);
            mat.data[i] = (Elem64)VSPIR_LHE_PRIME - magnitude;
        }
        //if ((Elem)mat.data[i] > error_norm) {
        //    error_norm = (Elem)mat.data[i];
        //}
    }
    //std::cout << "error norm = " << error_norm << std::endl;
}

template <typename MatrixType>
void constant(MatrixType& mat, const Elem val) {
    size_t len = mat.rows * mat.cols;
    for (size_t i = 0; i < len; i++) {
        mat.data[i] = val;
    }
}

template void constant(Matrix& mat, const Elem val);
template void constant(BinaryMatrix& mat, const Elem val);

bool eq(const Matrix& a, const Matrix& b, const bool verbose) {
    if (a.rows != b.rows || a.cols != b.cols) {
        if (verbose) {
            std::cout << "Dimension mismatch!\n";
            std::cout << a.rows << " x " << a.cols << " vs. " << b.rows << " x " << b.cols << std::endl;
        } 
        return false;
    }
    size_t len = a.rows * b.cols;
    for (size_t i = 0; i < len; i++) {
        if (a.data[i] != b.data[i]) {
            if (verbose) {
                std::cout << "Data mismatch!\n";
                std::cout << i << " " << a.data[i] << " " << b.data[i] << std::endl;
            }
            return false;
        }
    }
    return true;
}

template <typename MatType>
void print(const MatType& mat) {
    for (size_t r = 0; r < mat.rows; r++) {
        for (size_t c = 0; c < mat.cols; c++) {
            std::cout << mat.data[r*mat.cols + c] << " ";
        }
        std::cout << std::endl;
    }
}

template void print(const Matrix& mat);
template void print(const BinaryMatrix& mat); 

void matAddInPlace(Matrix& a, const Matrix& b, const Elem modulus) {
    const size_t rows = a.rows;
    const size_t cols = a.cols;

    if (rows != b.rows || cols != b.cols) {
        std::cout << "dimension mismatch\n!";
        assert(false);
    }

    const size_t len = rows * cols;
    ElemDouble big_temp;

    if (modulus == 0) {
        for (size_t i = 0; i < len; i++) {
            big_temp = ((ElemDouble)a.data[i]) + ((ElemDouble)b.data[i]);
            big_temp = fast_mod_max_64(static_cast<uint64_t>(big_temp));
            a.data[i] = static_cast<Elem>(big_temp);
            //a.data[i] += b.data[i];
        }
    } else {
        for (size_t i = 0; i < len; i++) {
            big_temp = ((ElemDouble)a.data[i]) + ((ElemDouble)b.data[i]);
            big_temp = fast_mod_max_64(static_cast<uint64_t>(big_temp));
            a.data[i] = static_cast<Elem>(big_temp);
            //a.data[i] += b.data[i];
            a.data[i] %= modulus;
        }
    }
}

void matAddInPlace64(Matrix64& a, const Matrix64& b, const Elem64 modulus) {
    const size_t rows = a.rows;
    const size_t cols = a.cols;

    if (rows != b.rows || cols != b.cols) {
        std::cout << "dimension mismatch\n!";
        assert(false);
    }

    const size_t len = rows * cols;
    absu128 big_temp;

    if (modulus == 0) {
        for (size_t i = 0; i < len; i++) {
            big_temp = ((absu128)a.data[i]) + ((absu128)b.data[i]);
            big_temp = fast_mod_vspir_max_64(static_cast<uint64_t>(big_temp));
            a.data[i] = static_cast<Elem64>(big_temp);
            //a.data[i] += b.data[i];
        }
    } else {
        for (size_t i = 0; i < len; i++) {
            big_temp = ((absu128)a.data[i]) + ((absu128)b.data[i]);
            big_temp = fast_mod_vspir_max_64(static_cast<uint64_t>(big_temp));
            a.data[i] = static_cast<Elem64>(big_temp);
            //a.data[i] += b.data[i];
            a.data[i] %= modulus;
        }
    }
}

Matrix matAdd(const Matrix& a, const Matrix& b) {
    const size_t rows = a.rows;
    const size_t cols = a.cols;

    if (rows != b.rows || cols != b.cols) {
        std::cout << "Input dimension mismatch\n!";
        assert(false);
    }

    Matrix out; out.init_no_memset(rows, cols);

    const size_t len = rows * cols;
    for (size_t i = 0; i < len; i++) {
        out.data[i] = a.data[i] + b.data[i];
    }

    return out;
}

Matrix matSub(const Matrix& a, const Matrix& b, const Elem modulus) {
    const size_t rows = a.rows;
    const size_t cols = a.cols;

    if (rows != b.rows || cols != b.cols) {
        std::cout << "Input dimension mismatch\n!";
        assert(false);
    }

    Matrix out; out.init_no_memset(rows, cols);

    const size_t len = rows * cols;
    ElemDouble big_temp;

    if (modulus == 0) {
        for (size_t i = 0; i < len; i++) {
            big_temp = ((ElemDouble)a.data[i]) + (ElemDouble)SOL_NTT_PRIME - ((ElemDouble)b.data[i]);
            big_temp = fast_mod_max_64(static_cast<uint64_t>(big_temp));
            out.data[i] = static_cast<Elem>(big_temp);
            //out.data[i] = a.data[i] - b.data[i];
        }
    } else {
        for (size_t i = 0; i < len; i++) {
            out.data[i] = a.data[i] + modulus - b.data[i];
            out.data[i] -= (out.data[i] >= modulus) ? modulus : 0;
        }
    }

    return out;
}

Matrix64 matSub64(const Matrix64& a, const Matrix64& b, const Elem64 modulus) {
    const size_t rows = a.rows;
    const size_t cols = a.cols;

    if (rows != b.rows || cols != b.cols) {
        std::cout << "Input dimension mismatch\n!";
        assert(false);
    }

    Matrix64 out; out.init_no_memset(rows, cols);

    const size_t len = rows * cols;
    absu128 big_temp;

    if (modulus == 0) {
        for (size_t i = 0; i < len; i++) {
            big_temp = ((absu128)a.data[i]) + (absu128)VSPIR_LHE_PRIME - ((absu128)b.data[i]);
            big_temp = fast_mod_vspir_max_64(static_cast<uint64_t>(big_temp));
            out.data[i] = static_cast<Elem64>(big_temp);
            //out.data[i] = a.data[i] - b.data[i];
        }
    } else {
        for (size_t i = 0; i < len; i++) {
            out.data[i] = a.data[i] + modulus - b.data[i];
            out.data[i] -= (out.data[i] >= modulus) ? modulus : 0;
        }
    }

    return out;
}

Matrix matMulScalar(const Matrix& a, const Elem b, const Elem modulus) {
    Matrix out; out.init_no_memset(a.rows, a.cols);
    const size_t len = a.rows * a.cols;
    if (modulus == 0) {
        for (size_t i = 0; i < len; i++) {
            out.data[i] = a.data[i] * b;
        }
    } else {
        for (size_t i = 0; i < len; i++) {
            out.data[i] = a.data[i] * b;
            out.data[i] %= modulus;
        }   
    }
    return out;
}

Matrix64 matMulScalar64(const Matrix64& a, const Elem64 b, const Elem64 modulus) {
    Matrix64 out; out.init_no_memset(a.rows, a.cols);
    const size_t len = a.rows * a.cols;
    if (modulus == 0) {
        for (size_t i = 0; i < len; i++) {
            out.data[i] = a.data[i] * b;
        }
    } else {
        for (size_t i = 0; i < len; i++) {
            out.data[i] = a.data[i] * b;
            out.data[i] %= modulus;
        }   
    }
    return out;
}

// Performs the divide and round step
Matrix matDivScalar(const Matrix& a, const Elem b) {
    Matrix out; out.init_no_memset(a.rows, a.cols);
    const size_t len = a.rows * a.cols;
    for (size_t i = 0; i < len; i++) {
        out.data[i] = a.data[i] / b;
        if (a.data[i] % b >= b/2) {
            out.data[i] += 1;
        }
    }
    return out;
}

// Performs the divide and round step
Matrix matDivScalar64(const Matrix64& a, const Elem64 b) {
    Matrix out; out.init_no_memset(a.rows, a.cols);
    const size_t len = a.rows * a.cols;
    for (size_t i = 0; i < len; i++) {
        Elem64 tmp = a.data[i] / b;
        out.data[i] = static_cast<Elem>(tmp);
        if (a.data[i] % b >= b/2) {
            out.data[i] += 1;
        }
    }
    return out;
}

// Performs multipy with c, then divide by b and round step
Matrix matMultThenDivScalar(const Matrix& a, const Elem b, const Elem c) {
    Matrix out; out.init_no_memset(a.rows, a.cols);
    const size_t len = a.rows * a.cols;
    ElemDouble big_temp, big_temp2;
    for (size_t i = 0; i < len; i++) {
        big_temp = ((ElemDouble)a.data[i]) * ((ElemDouble)c);
        big_temp2 = big_temp / ((ElemDouble)b);
        out.data[i] = static_cast<Elem>(big_temp2);
        //out.data[i] = a.data[i] / b;
        //if (a.data[i] % b >= b/2) {
        if (big_temp % ((ElemDouble)b) >= ((ElemDouble)b)/2) {
            out.data[i] += 1;
        }
        out.data[i] %= c;
    }
    return out;
}

Matrix matMul(const Matrix& a, const Matrix& b, const Elem modulus) {
    const size_t aRows = a.rows;
    const size_t aCols = a.cols;
    const size_t bCols = b.cols;

    if (aCols != b.rows) {
        std::cout << aCols << " " << b.rows << std::endl;
        std::cout << "Dimension mismatch!\n";
        assert(false);
    }

    //std::cout << "HERE NOW Rows: " << b.rows << " Cols: " << b.cols << std::endl;
    //for (uint64_t i = 0; i < b.rows; i++){
    //    for (uint64_t j = 0; j < b.cols; j++){
    //        if (j == 0) std::cout << " " << b.data[i*b.cols + j];
    //    }
    //}
    //std::cout << std::endl;

    Matrix out(aRows, bCols);  // memset values to zero
    ElemDouble big_temp;

    if (modulus == 0) {
        for (size_t i = 0; i < aRows; i++) {
            for (size_t k = 0; k < aCols; k++) {
                for (size_t j = 0; j < bCols; j++) {
                    big_temp = ((ElemDouble)a.data[aCols * i + k]) * ((ElemDouble)b.data[bCols * k + j]);
                    big_temp = fast_mod(big_temp);
                    big_temp = big_temp + (ElemDouble)out.data[bCols * i + j];
                    big_temp = fast_mod_max_64(static_cast<uint64_t>(big_temp));
                    out.data[bCols * i + j] = static_cast<Elem>(big_temp);
                    //out.data[bCols * i + j] += a.data[aCols * i + k] * b.data[bCols * k + j];
                    //if (i + j == 0) {
                    //    std::cout << "VSPIR a: " << a.data[aCols * i + k] << " b: " << b.data[bCols * k + j] << std::endl;
                    //}
                }
            }
        }

    } else {
        for (size_t i = 0; i < aRows; i++) {
            for (size_t k = 0; k < aCols; k++) {
                for (size_t j = 0; j < bCols; j++) {
                    out.data[bCols * i + j] += a.data[aCols * i + k] * b.data[bCols * k + j];
                    out.data[bCols * i + j] %= modulus;
                }
            }
        }

    }

    return out;
}

Matrix matMulSecondary(const Matrix& a, const Matrix& b, const Elem modulus) {
    const size_t aRows = a.rows;
    const size_t aCols = a.cols;
    const size_t bCols = b.cols;

    if (aCols != b.rows) {
        std::cout << aCols << " " << b.rows << std::endl;
        std::cout << "Dimension mismatch!\n";
        assert(false);
    }

    //std::cout << "HERE NOW Rows: " << b.rows << " Cols: " << b.cols << std::endl;
    //for (uint64_t i = 0; i < b.rows; i++){
    //    for (uint64_t j = 0; j < b.cols; j++){
    //        if (j == 0) std::cout << " " << b.data[i*b.cols + j];
    //    }
    //}
    //std::cout << std::endl;

    Matrix out(aRows, bCols);  // memset values to zero
    ElemDouble big_temp;

    if (modulus == 0) {
        for (size_t i = 0; i < aRows; i++) {
            for (size_t k = 0; k < aCols; k++) {
                for (size_t j = 0; j < bCols; j++) {
                    big_temp = ((ElemDouble)a.data[aCols * i + k]) * ((ElemDouble)b.data[bCols * k + j]);
                    big_temp = fast_mod_secondary(big_temp);
                    big_temp = big_temp + (ElemDouble)out.data[bCols * i + j];
                    big_temp = fast_mod_secondary_max_64(static_cast<uint64_t>(big_temp));
                    out.data[bCols * i + j] = static_cast<Elem>(big_temp);
                    //out.data[bCols * i + j] += a.data[aCols * i + k] * b.data[bCols * k + j];
                    //if (i + j == 0) {
                    //    std::cout << "VSPIR a: " << a.data[aCols * i + k] << " b: " << b.data[bCols * k + j] << std::endl;
                    //}
                }
            }
        }

    } else {
        for (size_t i = 0; i < aRows; i++) {
            for (size_t k = 0; k < aCols; k++) {
                for (size_t j = 0; j < bCols; j++) {
                    out.data[bCols * i + j] += a.data[aCols * i + k] * b.data[bCols * k + j];
                    out.data[bCols * i + j] %= modulus;
                }
            }
        }

    }

    return out;
}

// produces the product a*[b | c], where c is a vector
// assume that c is a vector
Matrix matMulAppendVec(const Matrix &a, const Matrix &b, const Matrix &c) {
    const size_t aRows = a.rows;
    const size_t aCols = a.cols;
    const size_t bCols = b.cols;

    if (aCols != b.rows) {
        std::cout << "Dimension mismatch!\n";
        std::cout << "trying to multiply ";
        std::cout << a.rows << " x " << a.cols << " times " << b.rows << " x " << b.cols << std::endl; 
        assert(false);
    }

    Matrix out(aRows, bCols + c.cols);  // memset values to zero

    for (size_t i = 0; i < aRows; i++) {
        for (size_t k = 0; k < aCols; k++) {
            const Elem val = a.data[aCols * i + k];
            for (size_t j = 0; j < bCols; j++) 
                out.data[out.cols * i + j] += val * b.data[bCols * k + j];
            for (size_t j = 0; j < c.cols; j++) 
                out.data[out.cols * i + bCols + j] += val * c.data[c.cols * k + j];
        }
    }

    return out;
}

Matrix binaryMatMulAppendVec(const BinaryMatrix &binary, const Matrix &b, const Matrix &c) {
    const size_t aRows = binary.rows;
    const size_t aCols = binary.cols;
    const size_t bCols = b.cols;

    if (aCols != b.rows) {
        std::cout << "Dimension mismatch!\n";
        assert(false);
    }

    Matrix out(aRows, bCols + c.cols);  // memset values to zero

    for (size_t i = 0; i < aRows; i++) {
        for (size_t k = 0; k < aCols; k++) {
            if (binary.data[aCols * i + k]) {
                for (size_t j = 0; j < bCols; j++) 
                    out.data[out.cols * i + j] += b.data[bCols * k + j];
                
                for (size_t j = 0; j < c.cols; j++) 
                    out.data[out.cols * i + bCols + j] += c.data[c.cols * k + j];
            }
        }
    }

    return out;
}

Matrix matMulLeftBinary(const BinaryMatrix& binary, const Matrix& b, const Elem modulus) {
    const size_t aRows = binary.rows;
    const size_t aCols = binary.cols;
    const size_t bCols = b.cols;

    if (aCols != b.rows) {
        std::cout << "Dimension mismatch!\n";
        assert(false);
    }

    Matrix out(aRows, bCols);  // memset values to zero
    ElemDouble big_temp;

    if (modulus == 0) {

        for (size_t i = 0; i < aRows; i++) {
            for (size_t k = 0; k < aCols; k++) {
                if (binary.data[aCols*i + k]) {
                    for (size_t j = 0; j < bCols; j++) {
                        big_temp = ((ElemDouble)out.data[bCols * i + j]) + ((ElemDouble)b.data[bCols * k + j]);
                        big_temp = fast_mod_max_64(static_cast<uint64_t>(big_temp));
                        out.data[bCols * i + j] = static_cast<Elem>(big_temp);
                        //out.data[bCols * i + j] += b.data[bCols * k + j];
                    }
                }
            }
        }

    } else {

        for (size_t i = 0; i < aRows; i++) {
            for (size_t k = 0; k < aCols; k++) {
                if (binary.data[aCols*i + k]) {
                    for (size_t j = 0; j < bCols; j++) {
                        out.data[bCols * i + j] += b.data[bCols * k + j];
                        out.data[bCols * i + j] %= modulus;
                    }
                }
            }
        }

    }

    return out;
}

Matrix matMulLeftBinarySecondary(const BinaryMatrix& binary, const Matrix& b, const Elem modulus) {
    const size_t aRows = binary.rows;
    const size_t aCols = binary.cols;
    const size_t bCols = b.cols;

    if (aCols != b.rows) {
        std::cout << "Dimension mismatch!\n";
        assert(false);
    }

    Matrix out(aRows, bCols);  // memset values to zero
    ElemDouble big_temp;

    if (modulus == 0) {

        for (size_t i = 0; i < aRows; i++) {
            for (size_t k = 0; k < aCols; k++) {
                if (binary.data[aCols*i + k]) {
                    for (size_t j = 0; j < bCols; j++) {
                        big_temp = ((ElemDouble)out.data[bCols * i + j]) + ((ElemDouble)b.data[bCols * k + j]);
                        big_temp = fast_mod_secondary_max_64(static_cast<uint64_t>(big_temp));
                        out.data[bCols * i + j] = static_cast<Elem>(big_temp);
                        //out.data[bCols * i + j] += b.data[bCols * k + j];
                    }
                }
            }
        }

    } else {

        for (size_t i = 0; i < aRows; i++) {
            for (size_t k = 0; k < aCols; k++) {
                if (binary.data[aCols*i + k]) {
                    for (size_t j = 0; j < bCols; j++) {
                        out.data[bCols * i + j] += b.data[bCols * k + j];
                        out.data[bCols * i + j] %= modulus;
                    }
                }
            }
        }

    }

    return out;
}

Matrix matMulRightBinary(const Matrix& b, const BinaryMatrix& binary, const Elem modulus) {
    Matrix rightMat = binary.asMatrix();
    return matMul(b, rightMat);
}

Matrix matMulRightBinarySecondary(const Matrix& b, const BinaryMatrix& binary, const Elem modulus) {
    Matrix rightMat = binary.asMatrix();
    return matMulSecondary(b, rightMat);
}

Matrix matMulVec(const Matrix& a, const Matrix& b, const Elem modulus) {
    const size_t aRows = a.rows;
    const size_t aCols = a.cols;
    const size_t bRows = b.rows;

    if (aCols != bRows || b.cols != 1) {
        std::cout << "Input dimension mismatch!\n";
        assert(false);
    }

    Matrix out; out.init_no_memset(aRows, 1);

    Elem tmp;
    ElemDouble big_tmp;

    if (modulus == 0) {  // machine word modulus

        /*std::cout << "Query u ";
        for (size_t j = 0; j < aCols; j++)
        {
            std::cout << b.data[j] << " ";
        }
        std::cout << std::endl;

        std::cout << "First row of DB ";
        for (size_t j = 0; j < aCols; j++)
        {
            std::cout << a.data[j] << " ";
        }
        std::cout << std::endl;*/

        for (size_t i = 0; i < aRows; i++)
        {
            tmp = 0;
            for (size_t j = 0; j < aCols; j++)
            {
                big_tmp = ((ElemDouble)a.data[aCols * i + j]) * ((ElemDouble)b.data[j]);
                big_tmp = fast_mod(big_tmp);
                big_tmp = big_tmp + (ElemDouble)tmp;
                big_tmp = fast_mod_max_64(static_cast<uint64_t>(big_tmp));
                tmp = static_cast<Elem>(big_tmp);
                //tmp += (Elem)big_tmp;
            }
            out.data[i] = tmp;
        }

    } else {

        for (size_t i = 0; i < aRows; i++)
        {
            tmp = 0;
            for (size_t j = 0; j < aCols; j++)
            {
                big_tmp = ((ElemDouble)a.data[aCols * i + j]) * ((ElemDouble)b.data[j]);
                big_tmp = fast_mod(big_tmp);
                big_tmp = big_tmp + (ElemDouble)tmp;
                big_tmp = fast_mod_max_64(static_cast<uint64_t>(big_tmp));
                tmp = static_cast<Elem>(big_tmp);
                //tmp += (Elem)big_tmp;
            }
            out.data[i] = tmp % modulus;
        }

    }

    return out;
}

Matrix64 matMulVec64(const Matrix64& a, const Matrix64& b, const Elem64 modulus) {
    const size_t aRows = a.rows;
    const size_t aCols = a.cols;
    const size_t bRows = b.rows;

    if (aCols != bRows || b.cols != 1) {
        std::cout << "Input dimension mismatch!\n";
        assert(false);
    }

    Matrix64 out; out.init_no_memset(aRows, 1);

    Elem64 tmp;
    absu128 big_tmp;

    if (modulus == 0) {  // machine word modulus

        /*std::cout << "Query u ";
        for (size_t j = 0; j < aCols; j++)
        {
            std::cout << b.data[j] << " ";
        }
        std::cout << std::endl;

        std::cout << "First row of DB ";
        for (size_t j = 0; j < aCols; j++)
        {
            std::cout << a.data[j] << " ";
        }
        std::cout << std::endl;*/

        for (size_t i = 0; i < aRows; i++)
        {
            tmp = 0;
            for (size_t j = 0; j < aCols; j++)
            {
                big_tmp = ((absu128)a.data[aCols * i + j]) * ((absu128)b.data[j]);
                big_tmp = fast_mod_vspir(big_tmp);
                big_tmp = big_tmp + (absu128)tmp;
                big_tmp = fast_mod_vspir_max_64(static_cast<uint64_t>(big_tmp));
                tmp = static_cast<Elem64>(big_tmp);
                //tmp += (Elem)big_tmp;
            }
            out.data[i] = tmp;
        }

    } else {

        for (size_t i = 0; i < aRows; i++)
        {
            tmp = 0;
            for (size_t j = 0; j < aCols; j++)
            {
                big_tmp = ((absu128)a.data[aCols * i + j]) * ((absu128)b.data[j]);
                big_tmp = fast_mod_vspir(big_tmp);
                big_tmp = big_tmp + (absu128)tmp;
                big_tmp = fast_mod_vspir_max_64(static_cast<uint64_t>(big_tmp));
                tmp = static_cast<Elem64>(big_tmp);
                //tmp += (Elem)big_tmp;
            }
            out.data[i] = tmp % modulus;
        }

    }

    return out;
}

Matrix matBinaryMulVec(const BinaryMatrix& a, const Matrix& b) {
    const size_t aRows = a.rows;
    const size_t aCols = a.cols;
    const size_t bRows = b.rows;

    if (aCols != bRows || b.cols != 1) {
        std::cout << "Input dimension mismatch!\n";
        assert(false);
    }

    Matrix out; out.init_no_memset(aRows, 1);

    /*std::cout << "Response v ";
    for (size_t j = 0; j < aCols; j++)
    {
        std::cout << b.data[j] << " ";
    }
    std::cout << std::endl;*/

    Elem tmp;
    Elem64 big_tmp;
    for (size_t i = 0; i < aRows; i++)
    {
        tmp = 0;
        for (size_t j = 0; j < aCols; j++)
        {
            if (a.data[aCols * i + j]) {
                big_tmp = (Elem64)tmp + (Elem64)b.data[j];
                big_tmp = fast_mod_max_64(big_tmp);
                tmp = static_cast<Elem>(big_tmp);
            }
            //tmp += (a.data[aCols * i + j]) ? b.data[j] : 0;
        }
        out.data[i] = tmp;
    }

    return out;
}

Matrix RowVecMulMat(const Matrix& v, const Matrix& m) {
    const size_t mRows = m.rows;
    const size_t mCols = m.cols;
    const size_t vCols = v.cols;

    if (vCols != mRows || v.rows != 1) {
        std::cout << "Input dimension mismatch!\n";
        assert(false);
    }

    Matrix out(1, mCols);

    Elem tmp;
    ElemDouble big_tmp;

    for (size_t c = 0; c < mCols; c++)
    {
        tmp = 0;
        for (size_t r = 0; r < mRows; r++)
        {
            big_tmp = ((ElemDouble)m.data[mCols * r + c]) * ((ElemDouble)v.data[r]);
            big_tmp = fast_mod(big_tmp);
            big_tmp = big_tmp + (ElemDouble)tmp;
            big_tmp = fast_mod_max_64(static_cast<uint64_t>(big_tmp));
            tmp = static_cast<Elem>(big_tmp);
            //tmp += (Elem)big_tmp;
        }
        out.data[c] = tmp;
    }

    return out;
}

Matrix RowVecMulBinaryMat(const Matrix& v, const BinaryMatrix& m) {
    const size_t mRows = m.rows;
    const size_t mCols = m.cols;
    const size_t vCols = v.cols;

    if (vCols != mRows || v.rows != 1) {
        std::cout << "Input dimension mismatch!\n";
        assert(false);
    }

    Matrix out(1, mCols);

    Elem tmp;
    Elem64 big_tmp;

    for (size_t c = 0; c < mCols; c++)
    {
        tmp = 0;
        for (size_t r = 0; r < mRows; r++)
        {
            if (m.data[mCols * r + c]) {
                big_tmp = (Elem64)tmp + (Elem64)v.data[r];
                big_tmp = fast_mod_max_64(big_tmp);
                tmp = static_cast<Elem>(big_tmp);
            }
        }
        out.data[c] = tmp;
    }
    return out;
}