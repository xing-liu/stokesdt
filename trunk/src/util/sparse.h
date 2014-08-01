/**
 * @file   sparse.h
 * @brief  Sparse matrix format definition
 */

#ifndef SPARSE_H_
#define SPARSE_H_


#include "common.h"


namespace stokesdt {

namespace detail {
    
/** 
 * @struct  SparseMatrix
 * @brief   Stores sparse matrix stored in blocked CSR format.
 */
typedef struct SparseMatrix {
    /// the dimension of the block 
    int sizeb;
    /// <code>sizeb^2</code>
    int sizeb2;
    /// the number of nonzero blocks
    int nnzb;
    /// the allocated storage capacity
    int maxnnzb;
    /// the number of blocks in row dimension
    int nrowb;
    /// the array of nonzero values
    double *val;
    /// the array of block column indices
    int *colbidx;
    /// the array of block row pointer
    int *rowbptr;
} SparseMatrix;


/// Creates a sparse matrix
bool CreateSparseMatrix(const int nrowb,
                        const int maxnnzb,
                        const int sizeb,
                        SparseMatrix **p_spmat);

/// Resizes the sparse matrix to contain <code>maxnnzb</code> nonzero blocks
bool ResizeSparseMatrix(const int maxnnzb,
                        SparseMatrix *spmat);

/// Destroys the sparse matrix
void DestroySparseMatrix(SparseMatrix *spmat);

/// Computes sparse matrix-vector multiplication
void SpMV3x3(const SparseMatrix *spmat,
             const int nrhs,
             const double alpha,
             const int ldx,
             const double *x,
             const double beta,
             const int ldy,
             double *y);

/// Transforms a spase matrix into dense format
void SparseToDense(const SparseMatrix *spmat, const int ldm, double *mat);

} // namespace detail

} // namespace stokesdt


#endif // SPARSE_H_
