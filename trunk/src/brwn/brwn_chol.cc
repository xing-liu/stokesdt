/**
 * @file   brwn_chol.cc
 * @brief  BrwnChol class implementation
 */

#include <mkl.h>
#include "brwn_chol.h"
#include "mob_matrix.h"
#include "log.h"
#include "profile.h"


namespace stokesdt {

BrwnChol::BrwnChol(const int dim) : dim_(dim)
{

}


BrwnChol::~BrwnChol()
{
    detail::AlignFree(cholmat_);
}


bool BrwnChol::Init()
{
    // check arguments
    if (dim_ <= 0) {
        LOG_ERROR("The specified dimension is less than or equal to 0.\n");
        return false;
    }

    LOG(3, "\n        Initializes BrwnChol\n");
    LOG(3, "        --------------------\n");
    LOG(3, "Dim = %d\n", dim_);
    
    // allocate buffer for the cholesky factor L
    ldm_ = detail::PadLen(dim_, sizeof(double));   
    cholmat_ = (double *)detail::AlignMalloc(sizeof(double) * dim_ * ldm_);
    if (NULL == cholmat_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * dim_ * ldm_);
        return false;
    }

    return true;
}


void BrwnChol::Compute(MobBase *mob, const int num_rhs,
                       const int ldz, const double *z,
                       const int ldy, double *y)
{
    START_TIMER(detail::BRWN_TICKS);
    
    // get the mobility matrix
    MobMatrix *mobmatrix = (MobMatrix *)mob;
    mobmatrix->GetMatrix(ldm_, cholmat_);

    // compute Cholesky factor L
    LAPACKE_dpotrf (LAPACK_COL_MAJOR, 'L', dim_, cholmat_, ldm_);

    // y = L' * z
    if (num_rhs == 1) {
        cblas_dcopy(dim_, z, 1, y, 1);
        cblas_dtrmv(CblasColMajor, CblasLower,
                    CblasNoTrans, CblasNonUnit,
                    dim_, cholmat_, ldm_, y, 1);
    } else {
        if (ldz == ldy) {
            cblas_dcopy(ldz * num_rhs, z, 1, y, 1);
        } else {
            for (int i = 0; i < num_rhs; i++) {
                cblas_dcopy(dim_, &z[i * ldz], 1, &y[i * ldy], 1);
            }
        }
        cblas_dtrmm(CblasColMajor, CblasLeft, CblasLower,
                    CblasNoTrans, CblasNonUnit,
                    dim_, num_rhs, 1.0, cholmat_, ldm_, y, ldy); 
    }

    STOP_TIMER(detail::BRWN_TICKS);
}

} // namespace stokesdt