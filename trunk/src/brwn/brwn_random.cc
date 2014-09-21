/**
 * @file   brwn_random.cc
 * @brief  BrwnRandom class implementation
 */

#include <mkl.h>
#include "brwn_random.h"
#include "mob_matrix.h"
#include "log.h"
#include "profile.h"


namespace stokesdt {

BrwnRandom::BrwnRandom(const int dim, const int num_vecs, const int max_nrhs)
    : dim_(dim), num_vecs_(num_vecs), max_nrhs_(max_nrhs)
{

}


BrwnRandom::~BrwnRandom()
{
    detail::AlignFree(u_);
    detail::AlignFree(s_);
    detail::AlignFree(v_);
    detail::AlignFree(tau_);
    detail::AlignFree(b_);
    detail::AlignFree(temp_);
    detail::AlignFree(first_term_);
    detail::AlignFree(second_term_);
    delete rnd_stream_;
}


bool BrwnRandom::Init()
{
    // check arguments
    if (dim_ <= 0) {
        LOG_ERROR("The specified dimension is less than or equal to 0.\n");
        return false;
    }
    if (num_vecs_ <= 0 || num_vecs_ > dim_) {
        LOG_ERROR("The specified number of random vectors"
                  " is less than or equal to 0 or larger than the"
                  " dimension of the mobility matrix.\n");
        return false;       
    }

    LOG(3, "\n        Initializes BrwnRandom\n");
    LOG(3, "        --------------------\n");
    LOG(3, "Dim      = %d\n", dim_);
    LOG(3, "Num-vecs = %d\n", num_vecs_);
    
    // allocate buffer for the random vectors
    ldm_ = detail::PadLen(dim_, sizeof(double));
    ldb_ = detail::PadLen(num_vecs_, sizeof(double));
    u_ = (double *)detail::AlignMalloc(sizeof(double) * num_vecs_ * ldm_);
    if (NULL == u_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * num_vecs_ * ldm_);
        return false;
    }
    s_ = (double *)detail::AlignMalloc(sizeof(double) * num_vecs_);
    if (NULL == s_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * num_vecs_);
        return false;
    }
    v_ = (double *)detail::AlignMalloc(sizeof(double) * num_vecs_ * ldm_);
    if (NULL == v_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * num_vecs_ * ldm_);
        return false;
    }
    tau_ = (double *)detail::AlignMalloc(sizeof(double) * num_vecs_);
    if (NULL == tau_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * num_vecs_);
        return false;
    }
    b_ = (double *)detail::AlignMalloc(sizeof(double) * num_vecs_ * ldb_);
    if (NULL == b_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * num_vecs_ * ldb_);
        return false;
    }
    temp_ = (double *)detail::AlignMalloc(sizeof(double) * max_nrhs_ * ldb_);
    if (NULL == temp_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * max_nrhs_ * ldb_);
        return false;
    }
    
    // allocate buffers for the first and second terms
    first_term_ =
        (double *)detail::AlignMalloc(sizeof(double) * max_nrhs_ * ldm_);
    if (NULL == first_term_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * max_nrhs_ * ldm_);
        return false;
    }
    second_term_ =
        (double *)detail::AlignMalloc(sizeof(double) * max_nrhs_ * ldm_);
    if (NULL == second_term_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * max_nrhs_ * ldm_);
        return false;
    }

    // init random number generator
    rnd_stream_ = new RndStream(SEED);
    if (!rnd_stream_->Init()) {
        return false;
    }
    
    return true;
}


void BrwnRandom::Compute(MobBase *mob, const int num_rhs,
                         const int ldz, const double *z,
                         const int ldy, double *y)
{
    if (mob != NULL) {        
        mob_ = mob;
        // compute omega
        rnd_stream_->Gaussian(0.0, 1.0, num_vecs_, ldm_, ldm_, u_);

        // v = M * omega
        mob_->MulVector(num_vecs_, 1.0, ldm_, u_, 0.0, ldm_, v_);
        mob_->MulVector(num_vecs_, 1.0, ldm_, v_, 0.0, ldm_, u_);
        mob_->MulVector(num_vecs_, 1.0, ldm_, u_, 0.0, ldm_, v_);
        
        // v = qr(v)        
        LAPACKE_dgeqrf(LAPACK_COL_MAJOR, dim_, num_vecs_, v_, ldm_, tau_);
        LAPACKE_dorgqr(LAPACK_COL_MAJOR, dim_, num_vecs_,
                       num_vecs_, v_, ldm_, tau_);
        // b = v' * M * v
        mob_->MulVector(num_vecs_, 1.0, ldm_, v_, 0.0, ldm_, u_);    
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                    num_vecs_, num_vecs_, dim_, 1.0, v_, ldm_, u_, ldm_,
                    0.0, b_, ldb_);
        // [uu ss ~] = svd(b), u = v*uu, s = sqrtm(ss)          
        LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'O', 'N',
                       num_vecs_, num_vecs_, b_, ldb_,
                       s_, NULL, num_vecs_, NULL, num_vecs_, tau_);
        for (int i = 0; i < num_vecs_; i++) {
            s_[i] = sqrt(s_[i]);
        }
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    dim_, num_vecs_, num_vecs_, 1.0, v_, ldm_, b_, ldb_,
                    0.0, u_, ldm_);
    }
    
    // compute brownian forces
    if (num_rhs == 1) {
        // temp = u'*z
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    num_vecs_, dim_, 1.0,
                    u_, ldm_, z, 1, 0.0, temp_, 1);
        // first_term = u*(s*temp)        
        for (int i = 0; i < num_vecs_; i++) {
            temp_[i] *= s_[i];
        }
        cblas_dgemv(CblasColMajor, CblasNoTrans,
                    dim_, num_vecs_, 1.0,
                    u_, ldm_, temp_, 1, 0.0, first_term_, 1);
        // second_term = z - v*(v'*z)
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    num_vecs_, dim_, 1.0,
                    v_, ldm_, z, 1, 0.0, temp_, 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans,
                    dim_, num_vecs_, 1.0,
                    v_, ldm_, temp_, 1, 0.0, second_term_, 1);
        cblas_daxpy(ldm_, -1.0, z, 1, second_term_, 1);       
        // beta2 = z'*M*z
        mob_->MulVector(1, 1.0, ldz, z, 0.0, ldy, y);
        double beta2 = cblas_ddot(dim_, z, 1, y, 1);
        // a = second_term'*second_term
        double a = cblas_ddot(dim_, second_term_, 1, second_term_, 1);
        // b = 2*first_term'*second_term
        double b = -2.0 * cblas_ddot(dim_, first_term_, 1, second_term_, 1);
        // c = first_term'*first_term - beta2
        double c = cblas_ddot(dim_, first_term_, 1, first_term_, 1) - beta2;
        // alpha = (-b + sqrt(b*b-4*a*c))/(2*a)
        double alpha = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
    
        // y = first_term + alpha*second_term
        cblas_dcopy(dim_, first_term_, 1, y, 1);
        cblas_daxpy(dim_, -alpha, second_term_, 1, y, 1);
    } else {
        // temp = u'*z
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    num_vecs_, num_rhs, dim_,
                    1.0, u_, ldm_, z, ldz,
                    0.0, temp_, ldb_);       
        // first_term = u*(s*temp)
        for (int i = 0; i < num_rhs; i++) {
            for (int j = 0; j < num_vecs_; j++) {
                temp_[i * ldb_ + j] *= s_[j];
            }
        }
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    dim_, num_rhs, num_vecs_,
                    1.0, u_, ldm_, temp_, ldb_,
                    0.0, first_term_, ldm_);
        // second_term = z - v*(v'*z)       
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    num_vecs_, num_rhs, dim_,
                    1.0, v_, ldm_, z, ldz,
                    0.0, temp_, ldb_);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    dim_, num_rhs, num_vecs_,
                    1.0, v_, ldm_, temp_, ldb_,
                    0.0, second_term_, ldm_);
        if (ldm_ == ldz) {        
            cblas_daxpy(ldm_ * num_rhs, -1.0, z, 1, second_term_, 1);
        } else {
            for (int k = 0; k < num_rhs; k++) {
                cblas_daxpy(ldm_, -1.0, &z[k * ldz], 1,
                            &second_term_[k * ldm_], 1);    
            }
        }

        mob_->MulVector(num_rhs, 1.0, ldz, z, 0.0, ldy, y);
        for (int k = 0; k < num_rhs; k++) {
            // beta2 = z'*M*z
            double beta2 = cblas_ddot(dim_, &z[k * ldz], 1, &y[k * ldy], 1);
            // a = second_term'*second_term
            double a = cblas_ddot(dim_, &second_term_[k * ldm_], 1,
                                  &second_term_[k * ldm_], 1);
            // b = 2*first_term'*second_term
            double b = -2.0 * cblas_ddot(dim_, &first_term_[k * ldm_],
                                         1, &second_term_[k * ldm_], 1);
            // c = first_term'*first_term - beta2
            double c = cblas_ddot(dim_, &first_term_[k * ldm_], 1,
                                  &first_term_[k * ldm_], 1) - beta2;
            // alpha = (-b + sqrt(b*b-4*a*c))/(2*a)
            double alpha = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
    
            // y = first_term + alpha*second_term
            cblas_dcopy(dim_, &first_term_[k * ldm_], 1, &y[k * ldy], 1);
            cblas_daxpy(dim_, -alpha, &second_term_[k * ldm_], 1,
                        &y[k * ldy], 1);
        } // for (int k = 0; k < num_rhs; k++)
    } // if (num_rhs == 1)
}

} // namespace stokesdt
