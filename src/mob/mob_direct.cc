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
#include "mob_direct.h"
#include "log.h"
#include "profile.h"


namespace stokesdt{

MobDirect::MobDirect(const int npos,
                     const double *rdi)
    : npos_(npos)
{

}


MobDirect::~MobDirect()
{
    detail::AlignFree(mat_);
}


bool MobDirect::Init()
{
    // check inputs
    if (npos_ <= 0) {
        LOG_ERROR("The specified number of particles"
                  " is less than or equal to 0: %d\n", npos_);
        return false;
    }

    LOG(3, "\n        Initializes MobDirect\n");
    LOG(3, "        ---------------------\n");
    
    // allocate buffer for the dense mobility matrix
    dim_mob_ = 3 * npos_;
    ldm_ = detail::PadLen(dim_mob_, sizeof(double));
    mat_ = (double *)detail::AlignMalloc(sizeof(double) * dim_mob_ * ldm_);
    if (NULL == mat_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n", 
            sizeof(double) * dim_mob_ * ldm_);
        return false;    
    }

    return true;
}


void MobDirect::Update(const double *pos, const double *rdi)
{
    START_TIMER(detail::MOB_TICKS);
    
    detail::NonEwaldKernel(npos_, pos, rdi, ldm_, mat_);
    
    STOP_TIMER(detail::MOB_TICKS);
}


void MobDirect::MulVector(const int num_rhs,
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


void MobDirect::GetMatrix(const int ldm, double *mat)
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