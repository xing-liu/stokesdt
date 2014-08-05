/**
 * @file   mob_ewald.cc
 * @brief  the MobEwald implementation
 */

#include <cstring>
#include <cstdio>
#include <stdint.h>
#include <vector>
#include <cmath>
#include <omp.h>
#include <mkl.h>
#include "mob_ewald.h"
#include "log.h"
#include "profile.h"


namespace stokesdt{

MobEwald::MobEwald(const int npos,
                   const double *rdi,
                   const double box_size,
                   const double tol,
                   const double xi)
    : npos_(npos),
      box_size_(box_size),
      tol_(tol),
      xi_(xi)
{

}



MobEwald::MobEwald(const int npos,
                   const double *rdi,
                   const double box_size,
                   const double tol)
    : npos_(npos),
      box_size_(box_size),
      tol_(tol)
{
    if (box_size > 0.0) {
        xi_ = pow(10.0, 1.0/6.0)*sqrt(M_PI) / box_size_;
    }
}


MobEwald::~MobEwald()
{
    detail::DestroyEwaldTable(ewald_tbl_);
    detail::AlignFree(mat_);
}


bool MobEwald::Init()
{
    // check inputs
    if (npos_ <= 0) {
        LOG_ERROR("The specified number of particles"
                  " is less than or equal to 0: %d\n", npos_);
        return false;
    }
    if (box_size_ <= 0.0) {
        LOG_ERROR("The specified dimension of the simulation box"
                  " is less than or equal to 0.0: %g\n", box_size_);
        return false;
    }
    if (tol_ <= 0.0) {
        LOG_ERROR("The requested tolerance"
                  " is less than or equal to 0.0: %g\n", tol_);
        return false;
    }
    if (xi_ <= 0.0) {
        LOG_ERROR("The specified Ewald paramter is less than"
                  " or equal to 0.0: %g\n", xi_);
        return false;
    }

    LOG(3, "\n        Initializes MobEwald\n");
    LOG(3, "        --------------------\n");
    LOG(3, "Box-size  = %g\n", box_size_);
    LOG(3, "Ewald-tol = %g\n", tol_);
    LOG(3, "Xi        = %g\n", xi_);
    
    // allocate buffer for the dense mobility matrix
    dim_mob_ = 3 * npos_;
    ldm_ = detail::PadLen(dim_mob_, sizeof(double));
    mat_ = (double *)detail::AlignMalloc(sizeof(double) * dim_mob_ * ldm_);
    if (NULL == mat_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n", 
            sizeof(double) * dim_mob_ * ldm_);
        return false;    
    }

    // precompute the ewald tables
    if (!detail::CreateEwaldTable(xi_, box_size_, tol_, &ewald_tbl_)) {
        return false;
    }

    return true;
}


void MobEwald::Update(const double *pos, const double *rdi)
{
    START_TIMER(detail::MOB_TICKS);
    
    if (ewald_tbl_ != NULL) {
        detail::EwaldKernel(xi_, ewald_tbl_, box_size_, npos_,
                            pos, rdi, ldm_, mat_, NULL, NULL);
    } else {
        LOG_WARN("The MobEwald is not properly initialized\n");
    }
    
    STOP_TIMER(detail::MOB_TICKS);
}


void MobEwald::MulVector(const int num_rhs,
                         const double alpha,
                         const int ldf,
                         const double *f,
                         const double beta,
                         const int ldv,
                         double *v)
{
    START_TIMER(detail::MOB_TICKS);
    
    if (num_rhs == 1){ // single vector
        cblas_dsymv(CblasRowMajor, CblasUpper,
                    dim_mob_, alpha, mat_, ldm_,
                    f, 1, beta, v, 1);
    } else{ // multiple vectors
        cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper,
                    dim_mob_, num_rhs, alpha, mat_, ldm_, f, ldf,
                    beta, v, ldv);            
    }

    STOP_TIMER(detail::MOB_TICKS);
}


void MobEwald::GetMatrix(const int ldm, double *mat)
{
    if (ldm == ldm_) {
        cblas_dcopy(dim_mob_ * ldm, mat_, 1, mat, 1);
    } else {
        for (int i = 0; i < dim_mob_; i++) {
            cblas_dcopy(dim_mob_, &(mat_[i * ldm_]), 1, &(mat[i * ldm]), 1);
        }
    }
}

} //namespace stokesdt