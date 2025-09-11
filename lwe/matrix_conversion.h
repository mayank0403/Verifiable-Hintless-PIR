#ifndef HINTLESS_PIR_MAT_CONV_H_
#define HINTLESS_PIR_MAT_CONV_H_

#include "types.h"
#include "verisimplepir/src/lib/pir/mat.h"
#include "verisimplepir/src/lib/pir/database.h"

// Name conflict resolution
typedef Matrix VMatrix;

// Converting between HintlessPIR and VeriSimplePIR matrix types.

void convertToVeriSimplePIRMatrix(VMatrix& vsp_matrix, const hintless_pir::lwe::Matrix& hintless_matrix);

void convertToVeriSimplePIRMatrix(VMatrix& vsp_matrix, const std::vector<std::vector<hintless_pir::lwe::Integer>>& hintless_matrix);

void convertToHintlessPIRMatrix(hintless_pir::lwe::Matrix& hintless_matrix, const VMatrix& vsp_matrix);

std::vector<std::vector<hintless_pir::lwe::Integer>> convertToHintlessPIRMatrix(const VMatrix& vsp_matrix, uint64_t rows, uint64_t cols);

std::vector<std::vector<hintless_pir::lwe::Integer>> convertToHintlessPIRMatrixColumnForm(const VMatrix& vsp_matrix, uint64_t rows, uint64_t cols);

std::vector<std::vector<hintless_pir::lwe::Integer64>> convertToHintlessPIRMatrixColumnForm64(const Matrix64& vsp_matrix, uint64_t rows, uint64_t cols);

std::vector<std::vector<hintless_pir::lwe::Integer>> convertToHintlessPIRMatricesFromVectors(const std::vector<VMatrix>& vsp_matrices, uint64_t rows);

std::vector<std::vector<hintless_pir::lwe::Integer64>> convertToHintlessPIRMatricesFromVectors64(const std::vector<Matrix64>& vsp_matrices, uint64_t rows);

void convertToVeriSimplePIRMatrixRowFormIn(VMatrix& vsp_matrix, const std::vector<std::vector<hintless_pir::lwe::Integer>>& hintless_matrix);

void convertToVeriSimplePIRMatrixRowFormIn64(Matrix64& vsp_matrix, const std::vector<std::vector<hintless_pir::lwe::Integer64>>& hintless_matrix);

std::vector<VMatrix> convertToVSPIRPIRVectors(const  std::vector<std::vector<hintless_pir::lwe::Integer>>& mat_in);

std::vector<Matrix64> convertToVSPIRPIRVectors64(const  std::vector<std::vector<hintless_pir::lwe::Integer64>>& mat_in);

VMatrix convertToVeriSimplePIRColumnMatrix(const hintless_pir::lwe::Vector& hintless_vec);

VMatrix convertToVeriSimplePIRColumnMatrixStdVec(const std::vector<hintless_pir::lwe::Integer>& hintless_vec);

VMatrix convertToVeriSimplePIRColumnMatrixSized(const hintless_pir::lwe::Vector& hintless_vec, uint64_t rows);

void convertToVectorFromVeriSimplePIRColumnMatrix(hintless_pir::lwe::Vector& hintless_vec, VMatrix& vsp_matrix);

void convertToVectorFromVeriSimplePIRColumnMatrix(std::vector<hintless_pir::lwe::Integer>& hintless_vec, VMatrix& vsp_matrix);

void convertToVeriSimplePIRDB(entry_t* vsp_db, const hintless_pir::lwe::Matrix& hintless_matrix, bool transpose = true);

void CompareMatrices(const VMatrix& vsp_matrix, const hintless_pir::lwe::Matrix& hintless_matrix, bool verbose);

void CompareMatrices(const VMatrix& vsp_matrix, const std::vector<std::vector<hintless_pir::lwe::Integer>>& hintless_matrix, bool verbose);

void CompareMatrices(const std::vector<VMatrix>& vsp_cols, const std::vector<std::vector<hintless_pir::lwe::Integer>>& hintless_matrix);

#endif  // HINTLESS_PIR_MAT_CONV_H_