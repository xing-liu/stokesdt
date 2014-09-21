/**
 * @file   brwn_lanczos.cc
 * @brief  BrwnLanczos class implementation
 */
 

#include <mkl.h>
#include <cstring>
#include "brwn_lanczos.h"
#include "log.h"
#include "profile.h"


namespace stokesdt {

int BrwnLanczos::Lanczos(MobBase *mob,
                         const double *z,
                         double *y)
{
    int ldh = detail::PadLen(max_iters_ + 1, sizeof(double));
    __declspec(align(detail::kAlignLen)) double d0[ldh];
    __declspec(align(detail::kAlignLen)) double e0[ldh];
    __declspec(align(detail::kAlignLen)) double d1[ldh];
    __declspec(align(detail::kAlignLen)) double e1[ldh];
    __declspec(align(detail::kAlignLen)) double h2[(max_iters_ + 1) * ldh];

    // v(:,1) = z/norm(z)
    double normz = cblas_dnrm2(dim_, z, 1);
    cblas_dcopy(dim_, z, 1, &(v_[0]), 1);
    cblas_dscal(dim_, 1.0/normz, &(v_[0]), 1);           
    int iter;
    for (iter = 0; iter < max_iters_; iter++) {
        // w = mat * v(i);  
        double *w = &(v_[(iter + 1) * ldv_]);      
        mob->MulVector(1, 1.0, ldv_, &(v_[iter * ldv_]), 0.0, ldv_, w);

        // w = w - h(i-1, i) * v(i-1)                
        if (iter > 0) {
            cblas_daxpy(dim_, -e0[iter - 1], &(v_[(iter - 1) * ldv_]), 1, w, 1);
        }

        // h(i, i) = <w, v[i]>
        d0[iter] = cblas_ddot(dim_, w, 1, &(v_[iter * ldv_]), 1);

        // w = w - h(i, i)*v(i)
        cblas_daxpy(dim_, -d0[iter], &(v_[iter * ldv_]), 1, w, 1);
        double normw = cblas_dnrm2(dim_, w, 1);
        e0[iter] = normw;

        // v(i + 1) = w = w/normw;
        cblas_dscal(dim_, 1.0/normw, w, 1);

        // h2 = sqrtm(h)      
        cblas_dcopy(iter + 1, d0, 1, d1, 1);
        cblas_dcopy(iter, e0, 1, e1, 1);
        LAPACKE_dstev(LAPACK_ROW_MAJOR, 'V', iter + 1, d1, e1, h2, ldh);
        #pragma vector aligned
        #pragma simd
        for (int j = 0; j < iter + 1; j++) {
            e1[j] = h2[j] * sqrt(d1[j]);
        }
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    iter + 1, iter + 1, 1.0, h2, ldh, e1, 1, 0.0, d1, 1);
            
        // y = v(:,1:i)*sqrth(:,1)*norm(z)
        cblas_dgemv(CblasColMajor, CblasNoTrans,
                    dim_, iter + 1, normz, v_, ldv_, d1, 1, 0.0, y, 1);

        // check convergence
        if (iter > 0) {
            cblas_daxpy(dim_, -1.0, y, 1, y_old_, 1);
            double normy = cblas_dnrm2(dim_, y_old_, 1);
            if (normy/normz < tol_) {
                break;
            }
        }
        cblas_dcopy(dim_, y, 1, y_old_, 1);
    } // for (iter = 0; i < max_iters_; i++)

    iter = (iter == max_iters_ ? -1 : iter);   
    return iter;
}


int BrwnLanczos::BlockLanczos(MobBase *mob,
                              const int num_rhs,
                              const int ldz,
                              const double *z,
                              const int ldy,
                              double *y)
{
    int ldr = detail::PadLen(num_rhs, sizeof(double));
    int ldh = detail::PadLen((max_iters_ + 1) * num_rhs, sizeof(double));    
    __declspec(align(detail::kAlignLen)) double band0[(num_rhs + 1) * ldh];
    __declspec(align(detail::kAlignLen)) double band1[(num_rhs + 1) * ldh];
    __declspec(align(detail::kAlignLen)) double D[(num_rhs + 1) * ldh];
    __declspec(align(detail::kAlignLen)) double H2[(max_iters_ + 1)*num_rhs*ldh];
    __declspec(align(detail::kAlignLen)) double H[num_rhs * ldr];
    __declspec(align(detail::kAlignLen)) double tau[num_rhs];
    __declspec(align(detail::kAlignLen)) double R[num_rhs * ldr];

    // [V(:,:,1), R] = qr(Z, 0)
    double normz = dlange("F", &dim_, &num_rhs, z, &ldz, NULL);
    for (int j = 0; j < num_rhs; j++) {
        cblas_dcopy(dim_, &(z[j * ldz]), 1, &(v_[j * ldv_]), 1);
    }
    LAPACKE_dgeqrf(LAPACK_COL_MAJOR, dim_, num_rhs, &(v_[0]), ldv_, tau);
    for (int j = 0; j < num_rhs; j++) {
        cblas_dcopy(num_rhs, &(v_[j * ldv_]), 1, &(R[j * ldr]), 1);
    }
    LAPACKE_dorgqr(LAPACK_COL_MAJOR, dim_, num_rhs, num_rhs, &v_[0], ldv_, tau);
    
    double beta = 0.0;
    int iter;
    for (iter = 0; iter < max_iters_; iter++) {
        // W = V(:,:,i-1) * H(:,:,i-1,i)
        double *W = &(v_[(iter + 1) * num_rhs * ldv_]);               
        if (iter > 0) {
            cblas_dcopy(num_rhs * ldv_,
                        &(v_[(iter - 1) * ldv_ * num_rhs]), 1, W, 1);
            cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper,
                        CblasTrans, CblasNonUnit,
                        dim_, num_rhs, 1.0, H, ldr, W, ldv_);
            beta = -1.0;
        } else {
            beta = 0.0;
        }
        // W = mat * V(:,:,i) - W;
        mob->MulVector(num_rhs, 1.0, ldv_, &(v_[iter * ldv_ * num_rhs]),
                       beta, ldv_, W);
        
        // H(:,:,i,i) = V(:,:,i)'*W;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    num_rhs, num_rhs, dim_,
                    1.0, &(v_[iter * ldv_ * num_rhs]), ldv_,
                    W, ldv_, 0.0, H, ldr);
     
        // W = W - V(:,:,i) * H(:,:,i,i)
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    dim_, num_rhs, num_rhs,
                    -1.0, &(v_[iter * ldv_ * num_rhs]), ldv_,
                    H, ldr, 1.0, W, ldv_);
                
        // [V(:,:,i+1),H(:,:,i+1,i)] = qr(W, 0)
        // H(:,:,i,i+1) = H(:,:,i+1,i)'
        LAPACKE_dgeqrf(LAPACK_COL_MAJOR, dim_, num_rhs, W, ldv_, tau);

        // H(:,:,i,i) = V(:,:,i)'*W;
        // H(:,:,i,i+1) = H(:,:,i+1,i)'
        for (int j = 0; j < num_rhs; j++) {
            for (int k = 0; k <= j; k++) {
                int row = j - k;
                int col = iter * num_rhs + k;
                band0[row * ldh + col] = H[j * ldr + k];                        
                int row2 = num_rhs - row;
                int col2 = iter * num_rhs + j;   
                band0[row2 * ldh + col2] = W[j * ldv_ + k];
            }
            cblas_dcopy(num_rhs, &(W[j * ldv_]), 1, &(H[j * ldr]), 1);
        }
        LAPACKE_dorgqr(LAPACK_COL_MAJOR, dim_, num_rhs, num_rhs, W, ldv_, tau);

        // H2 = H^.5;
        int kd = num_rhs;
        kd = (kd < (iter + 1) * num_rhs - 1 ? kd : (iter + 1) * num_rhs - 1);
        memcpy(band1, band0, sizeof(double) * ldh * (num_rhs + 1));
        LAPACKE_dsbev(LAPACK_ROW_MAJOR, 'V', 'L', (iter + 1) * num_rhs, kd,
                      band1, ldh, D, H2, ldh);
        for (int j = 0; j < num_rhs; j++) {
            #pragma simd
            for (int k = 0; k < (iter + 1) * num_rhs; k++) {
                band1[j * ldh + k] = H2[j * ldh + k] * sqrt(D[k]);
            }
        }               
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    (iter + 1)*num_rhs, num_rhs, (iter + 1)*num_rhs,
                    1.0, H2, ldh, band1, ldh, 0.0, D, ldh);
                 

        // Y = VV * H2(:, 1:k) * R;
        cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper,
                    CblasNoTrans, CblasNonUnit,
                    (iter + 1)*num_rhs, num_rhs, 1.0, R, ldr, D, ldh);    
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    dim_, num_rhs, num_rhs*(iter + 1), 1.0, v_, ldv_,
                    D, ldh, 0.0, y, ldv_);

        // check convergence
        if (iter > 0) {
            cblas_daxpy(ldv_ * num_rhs, -1.0, y, 1, y_old_, 1);
            double normy = dlange("F", &dim_, &num_rhs, y_old_, &ldv_, NULL);
            if (normy/normz < tol_) {
                break;
            }
        }
        cblas_dcopy(ldv_ * num_rhs, y, 1, y_old_, 1);
    } // for (iter = 0; i < max_iters_; i++)

    iter = (iter == max_iters_ ? -1 : iter); 
    return iter;
}


BrwnLanczos::BrwnLanczos(const int dim, const int max_iters,
                         const int max_nrhs, const double tol)
    : dim_(dim), max_iters_(max_iters), max_nrhs_(max_nrhs), tol_(tol)
{

}


BrwnLanczos::~BrwnLanczos()
{
    detail::AlignFree(v_);
    detail::AlignFree(y_old_); 
}


bool BrwnLanczos::Init()
{
    // check arguments
    if (dim_ <= 0) {
        LOG_ERROR("The specified dimension is less than or equal to 0.\n");
        return false;
    }
    if (max_iters_ <= 0) {
        LOG_ERROR("The specified maximum number of iterations"
                  " is less than or equal to 0.");
        return false;
    }
    if (max_nrhs_ <= 0) {
        LOG_ERROR("The specified maximum number of right hand sides"
                  " is less than or equal to 0.");
        return false;        
    }
    if (tol_ <= 0.0) {
        LOG_ERROR("The requested tolerance is less than or equal to 0.0.\n");
        return false;
    }
    if (max_nrhs_ > dim_/3) {
        LOG_WARN("The specified max_nrhs (%d) is too large."
                 " Set max_nrhs to %d.\n", max_nrhs_, dim_/3);        
        max_nrhs_ = dim_/3;       
    }

    LOG(3, "\n        Initializes BrwnLanczos\n");
    LOG(3, "        -----------------------\n");
    LOG(3, "Dim       = %d\n", dim_);
    LOG(3, "Max-iters = %d\n", max_iters_);
    LOG(3, "Max-nrhs  = %d\n", max_nrhs_);
    LOG(3, "Tol       = %g\n", tol_);

    // allocate buffer for the Krylov sequence
    ldv_ = detail::PadLen(dim_, sizeof(double));
    int len_v = max_nrhs_ * (max_iters_ + 1) * ldv_;
    v_ = (double *)detail::AlignMalloc(sizeof(double) * len_v);
    if (NULL == v_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n", sizeof(double) * len_v);
        return false;
    }

    y_old_ = (double *)detail::AlignMalloc(sizeof(double) * max_nrhs_ * ldv_);
    if (NULL == y_old_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * max_nrhs_ * ldv_);
        return false;
    }

    return true;
}


void BrwnLanczos::Compute(MobBase *mob, const int num_rhs,
                          const int ldz, const double *z,
                          const int ldy, double *y)
{
    START_TIMER(detail::BRWN_TICKS);
    
    // compute in groups
    for (int irhs = 0; irhs < num_rhs; irhs += max_nrhs_)
    {
        int nrhs = (irhs + max_nrhs_ > num_rhs ? num_rhs - irhs: max_nrhs_);
        int niters;        
        if (nrhs == 1) { // single Lanczos
            niters = Lanczos(mob, &(z[ldz * irhs]), &(y[ldy* irhs]));
        } else { // block Lanczos
            niters = BlockLanczos(mob, nrhs, ldz, &(z[ldz * irhs]),
                                  ldy, &(y[ldy* irhs]));
        }
        if (niters > 0) {
            LOG(3, "Computed Lanczos (%d vectors): converged in %d steps\n",
                nrhs, niters);
        } else {
            LOG(3, "Computed Lanczos (%d vectors): NOT converged\n", nrhs);
        }
    }

    START_TIMER(detail::BRWN_TICKS);
}

} // namespace stokesdt
