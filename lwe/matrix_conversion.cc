#include "matrix_conversion.h"
#include <iostream>

void convertToVeriSimplePIRMatrix(VMatrix& vsp_matrix, const hintless_pir::lwe::Matrix& hintless_matrix){
    // Assumes vsp_matrix is already allocated
    uint64_t rows, cols;
    rows = hintless_matrix.rows();
    cols = hintless_matrix.cols();
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            vsp_matrix.data[i*cols + j] = hintless_matrix(i, j);
        }
    }
    //std::cout << "Converted matrix from Hintless to VSP\n";
}

void convertToVeriSimplePIRMatrix(VMatrix& vsp_matrix, const std::vector<std::vector<hintless_pir::lwe::Integer>>& hintless_matrix){
    // Assumes vsp_matrix is already allocated; Hintless is organized by columns
    uint64_t rows, cols;
    cols = hintless_matrix.size();
    rows = hintless_matrix[0].size();
    //std::cout << "Rows: " << rows << " Cols: " << cols << std::endl;
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            vsp_matrix.data[i*cols + j] = hintless_matrix[j][i];
            //if (j == 0) std::cout << " " << hintless_matrix[j][i];
        }
    }
    //std::cout << std::endl;
    //std::cout << "Converted matrix from Hintless to VSP\n";
}

void convertToHintlessPIRMatrix(hintless_pir::lwe::Matrix& hintless_matrix, const VMatrix& vsp_matrix){
    // Assumes hintless_matrix is already allocated
    uint64_t rows, cols;
    rows = hintless_matrix.rows();
    cols = hintless_matrix.cols();
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            hintless_matrix(i, j) = vsp_matrix.data[i*cols + j];
        }
    }
    //std::cout << "Converted matrix from VSP to Hintless\n";
}

std::vector<std::vector<hintless_pir::lwe::Integer>> convertToHintlessPIRMatrix(const VMatrix& vsp_matrix, uint64_t rows, uint64_t cols){
    std::vector<std::vector<hintless_pir::lwe::Integer>> hintless(rows);
    for (uint64_t i = 0; i < rows; i++){
        hintless[i].resize(cols);
        for (uint64_t j = 0; j < cols; j++){
            hintless[i][j] = vsp_matrix.data[i*cols + j];
        }
    }
    //std::cout << "Converted matrix from VSP to Hintless\n";
    return hintless;
}

// Input and output in column format.
std::vector<std::vector<hintless_pir::lwe::Integer>> convertToHintlessPIRMatricesFromVectors(const std::vector<VMatrix>& vsp_matrices, uint64_t rows){
    uint32_t cols = vsp_matrices.size();
    std::vector<std::vector<hintless_pir::lwe::Integer>> hintless(cols);
    for (uint64_t i = 0; i < cols; i++){
        hintless[i].resize(rows);
        for (uint64_t j = 0; j < rows; j++){
            hintless[i][j] = vsp_matrices[i].data[j];
        }
    }
    //std::cout << "Converted matrix from VSP to Hintless\n";
    return hintless;
}

// Input and output in column format.
std::vector<std::vector<hintless_pir::lwe::Integer64>> convertToHintlessPIRMatricesFromVectors64(const std::vector<Matrix64>& vsp_matrices, uint64_t rows){
    uint32_t cols = vsp_matrices.size();
    std::vector<std::vector<hintless_pir::lwe::Integer64>> hintless(cols);
    for (uint64_t i = 0; i < cols; i++){
        hintless[i].resize(rows);
        for (uint64_t j = 0; j < rows; j++){
            hintless[i][j] = vsp_matrices[i].data[j];
        }
    }
    //std::cout << "Converted matrix from VSP to Hintless\n";
    return hintless;
}

// Input and output in column format.
std::vector<VMatrix> convertToVSPIRPIRVectors(const  std::vector<std::vector<hintless_pir::lwe::Integer>>& mat_in){
    uint32_t cols = mat_in.size();
    uint32_t rows = mat_in[0].size();
    std::vector<VMatrix> vsp_vectors(cols);
    
    for (uint64_t i = 0; i < cols; i++){
        vsp_vectors[i] = VMatrix(rows, 1);
        for (uint64_t j = 0; j < rows; j++){
            vsp_vectors[i].data[j] = mat_in[i][j];
        }
    }
    //std::cout << "Converted matrix to VSP from Hintless\n";
    return vsp_vectors;
}

// Input and output in column format.
std::vector<Matrix64> convertToVSPIRPIRVectors64(const  std::vector<std::vector<hintless_pir::lwe::Integer64>>& mat_in){
    uint32_t cols = mat_in.size();
    uint32_t rows = mat_in[0].size();
    std::vector<Matrix64> vsp_vectors(cols);
    
    for (uint64_t i = 0; i < cols; i++){
        vsp_vectors[i] = Matrix64(rows, 1);
        for (uint64_t j = 0; j < rows; j++){
            vsp_vectors[i].data[j] = mat_in[i][j];
        }
    }
    //std::cout << "Converted matrix to VSP from Hintless\n";
    return vsp_vectors;
}

std::vector<std::vector<hintless_pir::lwe::Integer>> convertToHintlessPIRMatrixColumnForm(const VMatrix& vsp_matrix, uint64_t rows, uint64_t cols){
    std::vector<std::vector<hintless_pir::lwe::Integer>> hintless(cols);
    for (uint64_t i = 0; i < cols; i++){
        hintless[i].resize(rows);
        for (uint64_t j = 0; j < rows; j++){
            hintless[i][j] = vsp_matrix.data[j*cols + i];
        }
    }
    //std::cout << "Converted matrix from VSP to Hintless\n";
    return hintless;
}

std::vector<std::vector<hintless_pir::lwe::Integer64>> convertToHintlessPIRMatrixColumnForm64(const Matrix64& vsp_matrix, uint64_t rows, uint64_t cols){
    std::vector<std::vector<hintless_pir::lwe::Integer64>> hintless(cols);
    for (uint64_t i = 0; i < cols; i++){
        hintless[i].resize(rows);
        for (uint64_t j = 0; j < rows; j++){
            hintless[i][j] = vsp_matrix.data[j*cols + i];
        }
    }
    //std::cout << "Converted matrix from VSP to Hintless\n";
    return hintless;
}

void convertToVeriSimplePIRMatrixRowFormIn(VMatrix& vsp_matrix, const std::vector<std::vector<hintless_pir::lwe::Integer>>& hintless_matrix){
    // Assumes vsp_matrix is already allocated; Hintless is organized by rows
    uint64_t rows, cols;
    rows = hintless_matrix.size();
    cols = hintless_matrix[0].size();
    //std::cout << "Rows: " << rows << " Cols: " << cols << std::endl;
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            vsp_matrix.data[i*cols + j] = hintless_matrix[i][j];
        }
    }
    //std::cout << std::endl;
    //std::cout << "Converted matrix from Hintless to VSP\n";
}

void convertToVeriSimplePIRMatrixRowFormIn64(Matrix64& vsp_matrix, const std::vector<std::vector<hintless_pir::lwe::Integer64>>& hintless_matrix){
    // Assumes vsp_matrix is already allocated; Hintless is organized by rows
    uint64_t rows, cols;
    rows = hintless_matrix.size();
    cols = hintless_matrix[0].size();
    //std::cout << "Rows: " << rows << " Cols: " << cols << std::endl;
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            vsp_matrix.data[i*cols + j] = hintless_matrix[i][j];
        }
    }
    //std::cout << std::endl;
    //std::cout << "Converted matrix from Hintless to VSP\n";
}

VMatrix convertToVeriSimplePIRColumnMatrix(const hintless_pir::lwe::Vector& hintless_vec){
    uint64_t rows, cols;
    rows = hintless_vec.rows();
    cols = 1;
    VMatrix vsp_matrix(rows, cols);
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            vsp_matrix.data[i*cols + j] = hintless_vec(i);
        }
    }
    //std::cout << "Converted vector from Hintless to VSP\n";
    return vsp_matrix;
}

VMatrix convertToVeriSimplePIRColumnMatrixStdVec(const std::vector<hintless_pir::lwe::Integer>& hintless_vec){
    uint64_t rows, cols;
    rows = hintless_vec.size();
    cols = 1;
    VMatrix vsp_matrix(rows, cols);
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            vsp_matrix.data[i*cols + j] = hintless_vec[i];
        }
    }
    //std::cout << "Converted vector from Hintless to VSP\n";
    return vsp_matrix;
}

VMatrix convertToVeriSimplePIRColumnMatrixSized(const hintless_pir::lwe::Vector& hintless_vec, uint64_t rows){
    uint64_t cols;
    assert(rows <= hintless_vec.rows());
    cols = 1;
    VMatrix vsp_matrix(rows, cols);
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            vsp_matrix.data[i*cols + j] = hintless_vec(i);
        }
    }
    //std::cout << "Converted vector from Hintless to VSP\n";
    return vsp_matrix;
}

void convertToVectorFromVeriSimplePIRColumnMatrix(hintless_pir::lwe::Vector& hintless_vec, VMatrix& vsp_matrix){
    // Assumes vsp_matrix is already allocated
    uint64_t rows, cols;
    rows = hintless_vec.rows();
    cols = 1;
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            hintless_vec(i) = vsp_matrix.data[i*cols + j];
        }
    }
    //std::cout << "Converted vector to Hintless from VSP\n";
}

void convertToVectorFromVeriSimplePIRColumnMatrix(std::vector<hintless_pir::lwe::Integer>& hintless_vec, VMatrix& vsp_matrix){
    // Assumes vsp_matrix is already allocated
    uint64_t rows, cols;
    rows = hintless_vec.size();
    cols = 1;
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            hintless_vec[i] = vsp_matrix.data[i*cols + j];
        }
    }
    //std::cout << "Converted vector to Hintless from VSP\n";
}

void convertToVeriSimplePIRDB(entry_t* vsp_db, const hintless_pir::lwe::Matrix& hintless_matrix, bool transpose){
    // Assumes vsp_db is already allocated
    // Transpose typically set to true because Hintless generates row-major, but VSP packs it in column-major
    uint64_t rows, cols;
    rows = hintless_matrix.rows();
    cols = hintless_matrix.cols();
    for (uint64_t i = 0; i < rows; i++){
        for (uint64_t j = 0; j < cols; j++){
            BigUnsigned val(0);
            val.setBlock(0, hintless_matrix(i, j));
            if (transpose) vsp_db[j*rows + i] = val;
            else
            vsp_db[i*cols + j] = val;
            if (i == 0 & j == 1) std::cout << " " << hintless_matrix(i, j) << " " << vsp_db[j*rows + i] << std::endl;
        }
    }
    //std::cout << "Converted matrix from Hintless to VSP\n";
}

void CompareMatrices(const VMatrix& vsp_matrix, const hintless_pir::lwe::Matrix& hintless_matrix, bool verbose){
    std::cout << "Comparing matrices init\n";
    uint64_t rows, cols;
    rows = hintless_matrix.rows();
    cols = hintless_matrix.cols();
    std::cout << "Comparing matrices\n";
    bool equal = true;
    uint64_t idx = -1;
    for (uint64_t i = 0; i < rows; i++){
        if (!equal) break;
        for (uint64_t j = 0; j < cols; j++){
            if(vsp_matrix.data[i*cols + j] != hintless_matrix(i, j)){
                equal = false;
                idx = i*cols + j;
                std::cout << vsp_matrix.data[i*cols + j] << " vs " << hintless_matrix(i, j) << std::endl;
                break;
            }
        }
    }
    if (equal){
        std::cout << "Matrices are equal\n";
    } else {
            std::cout << "Matrices are not equal at index " << idx << " out of " << rows*cols - 1 << std::endl;
            if (verbose){
            std::cout << "Left matrix:";
            for (uint64_t i = 0; i < rows; i++){
                for (uint64_t j = 0; j < cols; j++){
                    //std::cout << vsp_matrix.data[i*cols + j] << " ";
                }
                //std::cout << std::endl;
            }
            std::cout << "Right matrix:";
            for (uint64_t i = 0; i < rows; i++){
                for (uint64_t j = 0; j < cols; j++){
                    //std::cout << hintless_matrix(i, j) << " ";
                }
                //std::cout << std::endl;
            }
        }

        assert(false);
    }
}

void CompareMatrices(const VMatrix& vsp_matrix, const std::vector<std::vector<hintless_pir::lwe::Integer>>& hintless_matrix, bool verbose){
    uint64_t rows, cols;
    rows = hintless_matrix.size();
    cols = hintless_matrix[0].size();
    std::cout << "Rows: " << rows << " Cols: " << cols << std::endl;
    bool equal = true;
    uint64_t idx = -1;
    for (uint64_t i = 0; i < rows; i++){
        if (!equal) break;
        for (uint64_t j = 0; j < cols; j++){
            if(vsp_matrix.data[i*cols + j] != hintless_matrix[i][j]){
                equal = false;
                idx = i*cols + j;
                std::cout << vsp_matrix.data[i*cols + j] << " vs " << hintless_matrix[i][j] << std::endl;
                break;
            }
        }
    }
    if (equal){
        std::cout << "Matrices are equal\n";
    } else {
            std::cout << "Matrices are not equal at index " << idx << " out of " << rows*cols - 1 << std::endl;
            if (verbose){
            std::cout << "Left matrix:";
            for (uint64_t i = 0; i < rows; i++){
                for (uint64_t j = 0; j < cols; j++){
                    //std::cout << vsp_matrix.data[i*cols + j] << " ";
                }
                //std::cout << std::endl;
            }
            std::cout << "Right matrix:";
            for (uint64_t i = 0; i < rows; i++){
                for (uint64_t j = 0; j < cols; j++){
                    //std::cout << hintless_matrix[i][j] << " ";
                }
                //std::cout << std::endl;
            }
        }

        assert(false);
    }
}

// Both inputs are in column format
void CompareMatrices(const std::vector<VMatrix>& vsp_cols, const std::vector<std::vector<hintless_pir::lwe::Integer>>& hintless_matrix){
    uint64_t rows, cols;
    cols = hintless_matrix.size();
    rows = hintless_matrix[0].size();
    std::cout << "Rows: " << rows << " Cols: " << cols << std::endl;
    bool equal = true;
    uint64_t idx = -1;
    for (uint64_t i = 0; i < cols; i++){
        if (!equal) break;
        for (uint64_t j = 0; j < rows; j++){
            if(vsp_cols[i].data[j] != hintless_matrix[i][j]){
                equal = false;
                idx = i*rows + j;
                std::cout << vsp_cols[i].data[j] << " vs " << hintless_matrix[i][j] << std::endl;
                break;
            }
        }
    }
    if (equal){
        std::cout << "Matrices are equal\n";
    } else {
            std::cout << "Matrices are not equal at index " << idx << " out of " << rows*cols - 1 << std::endl;
        assert(false);
    }
}