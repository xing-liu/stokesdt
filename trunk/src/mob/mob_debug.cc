/**
 * @file   mob_debug.cc
 * @brief  MobDebug class implementation
 */

#include <cstring>
#include <cstdio>
#include <stdint.h>
#include <vector>
#include <cmath>
#include <omp.h>
#include <mkl.h>
#include "mob_debug.h"
#include "log.h"


namespace stokesdt{

MobDebug::MobDebug(const int npos,
                   const double *rdi,
                   const double box_size,
                   const double tol,
                   const double xi,
                   MobDebugType mode)
    : npos_(npos),
      box_size_(box_size),
      tol_(tol),
      xi_(xi),
      mode_(mode)
{

}


MobDebug::MobDebug(const int npos,
                   const double *rdi,
                   const double box_size,
                   const double tol,
                   MobDebugType mode)
    : npos_(npos),
      box_size_(box_size),
      tol_(tol),
      mode_(mode)
{
    if (box_size > 0.0) {
        xi_ = pow(10.0, 1.0/6.0)*sqrt(M_PI) / box_size_;
    }
}


MobDebug::~MobDebug()
{
    detail::DestroyEwaldTable(ewald_tbl_);
}


bool MobDebug::Init()
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
        LOG_ERROR("The specified requested tolerance"
                  " is less than or equal to 0.0: %g\n", tol_);
        return false;
    }
    if (xi_ <= 0.0) {
        LOG_ERROR("The specified Ewald paramter is less than"
                  " or equal to 0.0: %g\n", xi_);
        return false;
    }
                   
    dim_mob_ = 3 * npos_;
    ewald_tbl_ = new detail::EwaldTable;
    if (mode_ == EWALD) {
        if (!detail::InitRealTable(xi_, box_size_, tol_, ewald_tbl_) ||
            !detail::InitRecipTable(xi_, box_size_, tol_, ewald_tbl_)) {
            return false;
        }
    #if 0
        fprintf(stdout, "Ewald: rmax = %g, kmax = %g, nr = %d nk = %d\n",
            ewald_tbl_->rmax, ewald_tbl_->kmax,
            ewald_tbl_->nr, ewald_tbl_->nk);
    #endif
    } else if (mode_ == EWALD_REAL) {
        if (!detail::InitRealTable(xi_, box_size_, tol_, ewald_tbl_)) {
            return false;
        }
    #if 0    
        fprintf(stdout, "Ewald Real: rmax = %g, nr = %d\n",
            ewald_tbl_->rmax, ewald_tbl_->nr);
    #endif
    } else if (mode_ == EWALD_RECIP) {
        if (!detail::InitRecipTable(xi_, box_size_, tol_, ewald_tbl_)) {
            return false;
        }
    #if 0    
        fprintf(stdout, "Ewald Recip: kmax = %g, nk = %d\n",
            ewald_tbl_->kmax, ewald_tbl_->nk);
    #endif
    }
    
    return true;
}


void MobDebug::MulVector(const double *pos,
                         const double *rdi,
                         const int num_rhs,
                         const double alpha,
                         const int ldf,
                         const double *f,
                         const double beta,
                         const int ldv,
                         double *v)
{
    if (NULL == ewald_tbl_) {
        LOG_WARN("The MobDebug is not properly initialized\n");
        return;
    }
    
    if (mode_ == EWALD) {
        detail::EwaldVectorKernel(xi_, ewald_tbl_, box_size_, npos_,
                                  pos, rdi, alpha, ldf, f,
                                  beta, ldv, v, NULL, NULL);
    } else if (mode_ == EWALD_REAL) {
        detail::EwaldVectorKernel(xi_, ewald_tbl_, box_size_, npos_,
                                  pos, rdi, alpha, ldf, f,
                                  beta, ldv, NULL, v, NULL);   
    } else if (mode_ == EWALD_RECIP) {
        detail::EwaldVectorKernel(xi_, ewald_tbl_, box_size_, npos_,
                                  pos, rdi, alpha, ldf, f,
                                  beta, ldv, NULL, NULL, v);
    }
}


void MobDebug::GetMatrix(const double *pos, const double *rdi,
                         const int ldm, double *mat)
{
    if (NULL == ewald_tbl_) {
        LOG_WARN("The MobDebug is not properly initialized\n");
        return;
    }
    
    if (mode_ == EWALD) {
        detail::EwaldKernel(xi_, ewald_tbl_, box_size_, npos_,
                            pos, rdi, ldm, mat, NULL, NULL);
    } else if (mode_ == EWALD_REAL) {
        detail::EwaldKernel(xi_, ewald_tbl_, box_size_, npos_,
                            pos, rdi, ldm, NULL, mat, NULL); 
    } else if (mode_ == EWALD_RECIP) {
        detail::EwaldKernel(xi_, ewald_tbl_, box_size_, npos_,
                            pos, rdi, ldm, NULL, NULL, mat);
    }
}

} //namespace stokesdt