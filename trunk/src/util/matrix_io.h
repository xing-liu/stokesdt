/**
 * @file   matrix_io.h
 * @brief  the MatrixIO class definition
 */

#ifndef MATRIX_IO_H_
#define MATRIX_IO_H_


#include "common.h"


namespace stokesdt {

/** @class  MatrixIO
 *  @brief  Wraps matrix I/O operations
 */
class MatrixIO {
  public:
    /// Class constructor
    MatrixIO();
        
    /// Class deconstructor
    virtual ~MatrixIO();

    /** @brief  Reads a dense matrix from the specified text file
     *
     *  @param[in]  filename  the text file to be read
     *  @param[in]  nrows     the number of rows of the matrix
     *  @param[in]  ncols     the number of columns of the matrix
     *  @param[in]  ldm       the leading dimension of the matrix
     *  @param[out] mat       the pointer to the matrix
     */
    bool Read(const char* filename,
              const int nrows, const int ncols,
              const int ldm, double *mat);

    /** @brief  Writes the specified dense matrix into a text file
     *
     *  @param[in] nrows     the number of rows of the matrix
     *  @param[in] ncols     the number of columns of the matrix
     *  @param[in] ldm       the leading dimension of the matrix
     *  @param[in] mat       the pointer to the matrix
     *  @param[in] filename  the text file to be written
     */
    bool Write(const int nrows, const int ncols,
               const int ldm, const double *mat,
               const char* filename);

    /** @brief  Prints the specified dense matrix to stdout
     *
     *  @param[in] nrows  the number of rows of the matrix
     *  @param[in] ncols  the number of columns of the matrix
     *  @param[in] ldm    the leading dimension of the matrix
     *  @param[in] mat    the pointer to the matrix
     *  @param[in] name   the name of the matrix
     */
    bool Print(const int nrows, const int ncols,
               const int ldm, const double *mat,
               const char *name);
                     
    /** @brief  Calculates the Frobenius norm of
     *          the difference between two dense matrices
     *
     *  @param[in]  nrows  the number of rows of the matrix
     *  @param[in]  ncols  the number of columns of the matrix
     *  @param[in]  ldm1   the leading dimension of the matrix 1
     *  @param[in]  mat1   the pointer to the matrix 1
     *  @param[in]  ldm2   the leading dimension of the matrix 2
     *  @param[in]  mat2   the pointer to the matrix 2
     *  @param[out] error  the pointer to the return Frobenius norm
     */
    bool Compare(const int nrows, const int ncols,
                 const int ldm1, const double *mat1,
                 const int ldm2, const double *mat2,
                 double *error);
    
  private:
    DISALLOW_COPY_AND_ASSIGN(MatrixIO);
};

} // namespace stokesdt


#endif // MATRIX_IO_H_