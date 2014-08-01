/**
 * @file   matrix_io.cc
 * @brief  the MatrixIO class implementation
 */

#include <mkl.h>
#include <cstdio>
#include <cstring>
#include "matrix_io.h"
#include "log.h"


namespace stokesdt {

MatrixIO::MatrixIO()
{
}


MatrixIO::~MatrixIO()
{
}


bool MatrixIO::Write(const int nrows, const int ncols,
                     const int ldm, const double *mat,
                     const char* filename)
{
    FILE *fp = fopen(filename, "w+");
    if (fp == NULL) {
        LOG_ERROR("Open file %s failed\n", filename);
        return false;
    }

    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            fprintf(fp, "%le\n", mat[i * ldm + j]);
        }
    }

    fclose(fp);
    return true;
}


bool MatrixIO::Print(const int nrows, const int ncols,
                     const int ldm, const double *mat,
                     const char *name)
{
    fprintf(stdout, "Mat %s:\n", name);
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            fprintf(stdout, "%le ", mat[i * ldm + j]);
        }
        fprintf(stdout, "\n");
    }

    return true;
}


bool MatrixIO::Read(const char* filename,
                    const int nrows, const int ncols,
                    const int ldm, double *mat)
{
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        LOG_ERROR("Open file %s failed\n", filename);
        return false;
    }

    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            char str[1024];
            if(fgets(str, 1024, fp) != NULL) {
                sscanf(str, "%le", &mat[i * ldm + j]);
            } else {
                LOG_ERROR("Check the file size.\n");
                return false;
            }
        }
    }

    fclose(fp);

    return true;    
}


bool MatrixIO::Compare(const int nrows, const int ncols,
                       const int ldm1, const double *mat1,
                       const int ldm2, const double *mat2,
                       double *error)
{
    double *diff = (double *)malloc(nrows * ldm1 * sizeof (double));
    if (diff == NULL) {
        LOG_ERROR("Failed to allocate memory:\n");
        return false;
    }
    memcpy(diff, mat1, ldm1 * nrows * sizeof(double));
    if (ldm1 == ldm2) {
        cblas_daxpy(ldm1 * nrows, -1.0, mat2, 1, diff, 1);
    } else {
        for (int i = 0; i < nrows; i++) {
            cblas_daxpy(ncols, -1.0, &mat2[i * ldm2], 1, &diff[i * ldm1], 1);
        }
    }

    // dlange is col-major
    double norm0 = dlange("F", &ncols, &nrows, mat1, &ldm1, NULL);
    double norm_diff = dlange("F", &ncols, &nrows, diff, &ldm1, NULL);
    *error = norm_diff/norm0;

    free(diff);
    return true;
}

} // namespace stokesdt