/**
 * @file   ewald.cc
 * @brief  Ewald implementation
 */

#include <cstring>
#include <vector>
#include <cmath>
#include <omp.h>
#include "common.h"
#include "ewald.h"
#include "log.h"


namespace stokesdt {

namespace detail {

static double EwaldReal(const double r, const double xi)
{
    double xir = xi * r;
    double xir2 = xir * xir;
    double r2 = r * r;
    double v = erfc(xir) * (0.75+1.5/r2) / r +
           exp(-xir2) * xi / sqrt(M_PI) *
           (4.5 + 3.0*xir2 +
           (3.0 + xir2 * (14.0 + xir2 * (4.0 + xir2))) / r2);
    return v;
}


static double EwaldRecip(const double k, const double xi)
{
    double kx = k / xi;
    double kx2 = kx * kx;
    double k2 = k * k;
    double v = 6.0 * M_PI * exp(-kx2/4.0) *
        (1.0 + k2/3.0) / k2 * (1.0 + kx2 * (1.0/4.0 + kx2/8.0));
    return v;
}


static double XitoRmax(const double xi, const double tol,
                       const double box_size, const double safe_factor)
{
    double red_factor = 0.99;
    double dr = box_size;
   
    // rmax
    double xi0 = xi / safe_factor;
    double rmax = 10.0 * sqrt(-log(tol)) / xi;
    while (1){
        double ewald0 = EwaldReal(rmax * red_factor, xi0);
        double ewald1 = 6.0 * EwaldReal (rmax * red_factor + dr, xi0);
        if (ewald0 + ewald1 >= tol){
            break;
        }
        rmax *= red_factor;
    }

    return rmax;
}


static double XitoKmax(const double xi, const double tol,
                       const double box_size, const double safe_factor)
{
    double red_factor = 0.99;
    double dk = 2.0 * M_PI / box_size;

    double xi0 = xi * safe_factor;
    double kmax = 10.0 * sqrt(-log(tol)) * 2.0 * xi;
    while (1){
        double ewald0 = EwaldRecip(kmax * red_factor, xi0);
        double ewald1 = 6.0 * EwaldRecip (kmax * red_factor + dk, xi0);
        if (ewald0 + ewald1 >= tol){
            break;
        }
        kmax *= red_factor;
    }

    return kmax;
}


void EwaldKernel(const double xi, const EwaldTable *ewald_tbl,
                 const double box_size, const int npos,
                 const double *pos, const double *rdi,
                 const int ldmat, double *mat,
                 double *mat_real, double *mat_recip)
{    
    // compute tile size
    int align_dp = detail::kAlignLen/sizeof(double);
    int tile_size;
    int num_tiles;   
    for (int i = 1; i <= align_dp; i++) {
        if (0 == (3 * i) % align_dp) {
            tile_size = i;
        }
    }
    num_tiles = (npos + tile_size - 1)/tile_size;    
    int64_t num_tiles2 = (int64_t)num_tiles * num_tiles;
    
    // using OpenMP
    double xi2 = xi * xi;    
    double rmax2 = ewald_tbl->rmax * ewald_tbl->rmax;
    double xiaspi = xi / sqrt(M_PI);
    int nr = ewald_tbl->nr;
    int nk = ewald_tbl->nk;
    double *rlx_tbl = ewald_tbl->rlx_tbl;
    double *rly_tbl = ewald_tbl->rly_tbl;
    double *rlz_tbl = ewald_tbl->rlz_tbl;    
    double *ex_tbl = ewald_tbl->ex_tbl;
    double *ey_tbl = ewald_tbl->ey_tbl;
    double *ez_tbl = ewald_tbl->ez_tbl;
    double *k1_tbl = ewald_tbl->k1_tbl;
    double *k2_tbl = ewald_tbl->k2_tbl;
    double *k3_tbl = ewald_tbl->k3_tbl;
    double *k_tbl  = ewald_tbl->k_tbl;
    double *kexp_tbl = ewald_tbl->kexp_tbl;
    #pragma omp parallel
    {
        __declspec(align(detail::kAlignLen)) int lm0[nr + detail::kSimdWidth];        
        #pragma omp for schedule(dynamic)
        for (int64_t tile = 0; tile < num_tiles2; tile++) {
            int starti = (int)(tile/num_tiles) * tile_size;
            int endi = starti + tile_size;
            endi = endi < npos ? endi : npos;
            int startj = (int)(tile%num_tiles) * tile_size;
            int endj = startj + tile_size;
            endj = endj < npos ? endj : npos;
            if (starti > startj) {
                continue;
            }           
            for (int i = starti; i < endi; i++) {
                int i0 = i * 3;
                const double *pos1 = &pos[3 * i];
                double aa = rdi[i];
                double aa2 = aa * aa;
                for (int j = startj; j < endj; j++) {
                    if (i > j) {
                        continue;
                    }
                    int j0 = j * 3;
                    const double *pos2 = &pos[3 * j];
                    double ab = rdi[j];
                    double ab2 = ab * ab;
                    double xx0 = drem(pos2[0] - pos1[0], box_size);
                    double yy0 = drem(pos2[1] - pos1[1], box_size);
                    double zz0 = drem(pos2[2] - pos1[2], box_size);
                    // self part
                    double real0 = 0.0;
                    double real1 = 0.0;
                    double real2 = 0.0;
                    double real3 = 0.0;
                    double real4 = 0.0;
                    double real5 = 0.0;
                    if (mat != NULL || mat_real != NULL) {
                        // self
                        if (i == j) {
                            double self_a = 1.0 / aa - xiaspi * 
                                            (6.0 - 40.0 / 3.0 * xi2 * aa2);
                            real0 = self_a;
                            real3 = self_a;
                            real5 = self_a;                            
                        }
                        // real part
                        int nr0 = 0;
                        for (int m = 0; m < nr; m++) {
                            double xx = xx0 + rlx_tbl[m];
                            double yy = yy0 + rly_tbl[m];
                            double zz = zz0 + rlz_tbl[m];
                            double rr = xx * xx + yy * yy + zz * zz;
                            if (rr != 0.0 && rr <= rmax2) {
                                double r = sqrt(rr);
                                lm0[nr0++] = m;
                                if (r < (aa + ab)) {
                                    double xa = 
                                        -(1.5 - (aa2 + ab2) * 0.5 / rr)/r +
                                        2.0 / (aa + ab) -
                                        3.0 * r / 8.0 / (aa2 + ab2);
                                    double ya = 
                                        -(0.75 + (aa2 + ab2) * 0.25/rr)/r +
                                        2.0 / (aa + ab) -
                                        9.0 * r /16.0/(aa2 + ab2);
                                    double ex = xx / r;
                                    double ey = yy / r;
                                    double ez = zz / r;
                                    real0 += (xa - ya) * ex * ex + ya;
                                    real1 += (xa - ya) * ex * ey;
                                    real2 += (xa - ya) * ex * ez;
                                    real3 += (xa - ya) * ey * ey + ya;
                                    real4 += (xa - ya) * ey * ez;
                                    real5 += (xa - ya) * ez * ez + ya;
                                }
                            }
                        }
                        #pragma vector aligned
                        #pragma simd
                        for (int m = 0; m < nr0; m++) {
                            int mmm = lm0[m];
                            double xx = xx0 + rlx_tbl[mmm];
                            double yy = yy0 + rly_tbl[mmm];
                            double zz = zz0 + rlz_tbl[mmm];
                            double rr = xx * xx + yy * yy + zz * zz;
                            double r = sqrt(rr);
                            double ex = xx / r;
                            double ey = yy / r;
                            double ez = zz / r;
                            double xa;
                            double ya;
                            scalars_ewald_F(xi, r, aa, ab, &xa, &ya);
                            real0 += (xa - ya) * ex * ex + ya;
                            real1 += (xa - ya) * ex * ey;
                            real2 += (xa - ya) * ex * ez;
                            real3 += (xa - ya) * ey * ey + ya;
                            real4 += (xa - ya) * ey * ez;
                            real5 += (xa - ya) * ez * ez + ya;
                        }
                        if (mat_real != NULL) {
                            mat_real[(i0 + 0) * ldmat + j0 + 0] = real0;
                            mat_real[(i0 + 0) * ldmat + j0 + 1] = real1;
                            mat_real[(i0 + 0) * ldmat + j0 + 2] = real2;
                            mat_real[(i0 + 1) * ldmat + j0 + 0] = real1;
                            mat_real[(i0 + 1) * ldmat + j0 + 1] = real3;
                            mat_real[(i0 + 1) * ldmat + j0 + 2] = real4;
                            mat_real[(i0 + 2) * ldmat + j0 + 0] = real2;
                            mat_real[(i0 + 2) * ldmat + j0 + 1] = real4;
                            mat_real[(i0 + 2) * ldmat + j0 + 2] = real5;
                            if (j0 != i0) {
                                mat_real[(j0 + 0) * ldmat + i0 + 0] = real0;
                                mat_real[(j0 + 0) * ldmat + i0 + 1] = real1;
                                mat_real[(j0 + 0) * ldmat + i0 + 2] = real2;
                                mat_real[(j0 + 1) * ldmat + i0 + 0] = real1;
                                mat_real[(j0 + 1) * ldmat + i0 + 1] = real3;
                                mat_real[(j0 + 1) * ldmat + i0 + 2] = real4;
                                mat_real[(j0 + 2) * ldmat + i0 + 0] = real2;
                                mat_real[(j0 + 2) * ldmat + i0 + 1] = real4;
                                mat_real[(j0 + 2) * ldmat + i0 + 2] = real5;
                            }
                        }
                    } // end of if (mat != NULL && mat_real != NULL)
                    
                    // recip part
                    double recip0 = 0.0;
                    double recip1 = 0.0;
                    double recip2 = 0.0;
                    double recip3 = 0.0;
                    double recip4 = 0.0;
                    double recip5 = 0.0;
                    if (mat != NULL || mat_recip != NULL) {
                        #pragma vector aligned
                        #pragma simd
                        for (int m = 0; m < nk; m++) {
                            double ex   = ex_tbl[m];
                            double ey   = ey_tbl[m];
                            double ez   = ez_tbl[m];
                            double k1   = k1_tbl[m];
                            double k2   = k2_tbl[m];
                            double k3   = k3_tbl[m];
                            double k    = k_tbl[m];
                            double kexp = kexp_tbl[m];
                            double kk = k * k;
                            double ya = 
                                6.0 * (1.0 - kk * (aa2 + ab2) / 6.0) * kexp;
                            double cf = cos(k1 * xx0 + k2 * yy0 + k3 * zz0);
                            recip0 +=  cf * ya * (1.0 - ex * ex);
                            recip1 += -cf * ya * ex * ey;
                            recip2 += -cf * ya * ex * ez;
                            recip3 +=  cf * ya * (1.0 - ey * ey);
                            recip4 += -cf * ya * ey * ez;
                            recip5 +=  cf * ya * (1.0 - ez * ez);
                        }
                        if (mat_recip != NULL) {
                            mat_recip[(i0 + 0) * ldmat + j0 + 0] = recip0;
                            mat_recip[(i0 + 0) * ldmat + j0 + 1] = recip1;
                            mat_recip[(i0 + 0) * ldmat + j0 + 2] = recip2;
                            mat_recip[(i0 + 1) * ldmat + j0 + 0] = recip1;
                            mat_recip[(i0 + 1) * ldmat + j0 + 1] = recip3;
                            mat_recip[(i0 + 1) * ldmat + j0 + 2] = recip4;
                            mat_recip[(i0 + 2) * ldmat + j0 + 0] = recip2;
                            mat_recip[(i0 + 2) * ldmat + j0 + 1] = recip4;
                            mat_recip[(i0 + 2) * ldmat + j0 + 2] = recip5;
                            if (j0 != i0) {
                                mat_recip[(j0 + 0) * ldmat + i0 + 0] = recip0;
                                mat_recip[(j0 + 0) * ldmat + i0 + 1] = recip1;
                                mat_recip[(j0 + 0) * ldmat + i0 + 2] = recip2;
                                mat_recip[(j0 + 1) * ldmat + i0 + 0] = recip1;
                                mat_recip[(j0 + 1) * ldmat + i0 + 1] = recip3;
                                mat_recip[(j0 + 1) * ldmat + i0 + 2] = recip4;
                                mat_recip[(j0 + 2) * ldmat + i0 + 0] = recip2;
                                mat_recip[(j0 + 2) * ldmat + i0 + 1] = recip4;
                                mat_recip[(j0 + 2) * ldmat + i0 + 2] = recip5;
                            }
                        }
                    } // end of if (mat != NULL && mat_recip != NULL)
                    
                    if (mat != NULL) {
                        double ewald0 = real0 + recip0;
                        double ewald1 = real1 + recip1;
                        double ewald2 = real2 + recip2;
                        double ewald3 = real3 + recip3;
                        double ewald4 = real4 + recip4;
                        double ewald5 = real5 + recip5;
                        mat[(i0 + 0) * ldmat + j0 + 0] = ewald0;
                        mat[(i0 + 0) * ldmat + j0 + 1] = ewald1;
                        mat[(i0 + 0) * ldmat + j0 + 2] = ewald2;
                        mat[(i0 + 1) * ldmat + j0 + 0] = ewald1;
                        mat[(i0 + 1) * ldmat + j0 + 1] = ewald3;
                        mat[(i0 + 1) * ldmat + j0 + 2] = ewald4;
                        mat[(i0 + 2) * ldmat + j0 + 0] = ewald2;
                        mat[(i0 + 2) * ldmat + j0 + 1] = ewald4;
                        mat[(i0 + 2) * ldmat + j0 + 2] = ewald5;
                        if (j0 != i0) {
                            mat[(j0 + 0) * ldmat + i0 + 0] = ewald0;
                            mat[(j0 + 0) * ldmat + i0 + 1] = ewald1;
                            mat[(j0 + 0) * ldmat + i0 + 2] = ewald2;
                            mat[(j0 + 1) * ldmat + i0 + 0] = ewald1;
                            mat[(j0 + 1) * ldmat + i0 + 1] = ewald3;
                            mat[(j0 + 1) * ldmat + i0 + 2] = ewald4;
                            mat[(j0 + 2) * ldmat + i0 + 0] = ewald2;
                            mat[(j0 + 2) * ldmat + i0 + 1] = ewald4;
                            mat[(j0 + 2) * ldmat + i0 + 2] = ewald5;
                        }
                    } // if (mat != NULL)
                } // for (int j = startj; j < endj; j++)
            } // for (int i = starti; i < endi; i++)
        } // end of #pragma omp for schedule(dynamic)
    } // end of #pragma omp parallel
}


void EwaldVectorKernel(const double xi, const EwaldTable *ewald_tbl,
                       const double box_size, const int npos,
                       const double *pos, const double *rdi,
                       const double alpha, const int ldf, const double *f,
                       const double beta, const int ldv, double *v,
                       double *v_real, double *v_recip)
{
    // Using OpenMP
    double xi2 = xi * xi;    
    double rmax2 = ewald_tbl->rmax * ewald_tbl->rmax;
    double xiaspi = xi / sqrt(M_PI);
    int nr = ewald_tbl->nr;
    int nk = ewald_tbl->nk;
    double *rlx_tbl = ewald_tbl->rlx_tbl;
    double *rly_tbl = ewald_tbl->rly_tbl;
    double *rlz_tbl = ewald_tbl->rlz_tbl;    
    double *ex_tbl = ewald_tbl->ex_tbl;
    double *ey_tbl = ewald_tbl->ey_tbl;
    double *ez_tbl = ewald_tbl->ez_tbl;
    double *k1_tbl = ewald_tbl->k1_tbl;
    double *k2_tbl = ewald_tbl->k2_tbl;
    double *k3_tbl = ewald_tbl->k3_tbl;
    double *k_tbl  = ewald_tbl->k_tbl;
    double *kexp_tbl = ewald_tbl->kexp_tbl;
    if (rlx_tbl == NULL) {
        nr = 0;
    }
    #pragma omp parallel
    {
        __declspec(align(detail::kAlignLen)) int lm0[nr + detail::kSimdWidth];        
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < npos; i++) {
            int i0 = i * 3;
            const double *pos1 = &pos[3 * i];
            double aa = rdi[i];
            double aa2 = aa * aa;
            // initialization
            if (v_real != NULL) {
                v_real[i0 + 0] *= beta;
                v_real[i0 + 1] *= beta;
                v_real[i0 + 2] *= beta;
            }
            if (v_recip != NULL) {
                v_recip[i0 + 0] *= beta;
                v_recip[i0 + 1] *= beta;
                v_recip[i0 + 2] *= beta;
            }
            if (v != NULL) {
                v[i0 + 0] *= beta;
                v[i0 + 1] *= beta;
                v[i0 + 2] *= beta;
            }
            for (int j = 0; j < npos; j++) {
                int j0 = j * 3;
                const double *pos2 = &pos[3 * j];
                double ab = rdi[j];
                double ab2 = ab * ab;
                double xx0 = drem(pos2[0] - pos1[0], box_size);
                double yy0 = drem(pos2[1] - pos1[1], box_size);
                double zz0 = drem(pos2[2] - pos1[2], box_size);
                // self part
                double real0 = 0.0;
                double real1 = 0.0;
                double real2 = 0.0;
                double real3 = 0.0;
                double real4 = 0.0;
                double real5 = 0.0;
                if (v != NULL || v_real != NULL) {
                    // self
                    if (i == j) {
                        double self_a =
                            1.0 / aa - xiaspi * (6.0 - 40.0 / 3.0 * xi2 * aa2);
                        real0 = self_a;
                        real3 = self_a;
                        real5 = self_a;
                    }
                    
                    // real part
                    int nr0 = 0;
                    for (int m = 0; m < nr; m++) {
                        double xx = xx0 + rlx_tbl[m];
                        double yy = yy0 + rly_tbl[m];
                        double zz = zz0 + rlz_tbl[m];
                        double rr = xx * xx + yy * yy + zz * zz;
                        if (rr != 0.0 && rr <= rmax2) {
                            double r = sqrt(rr);
                            lm0[nr0++] = m;
                            if (r < (aa + ab)) {
                                double xa = -(1.5 - (aa2 + ab2) * 0.5 / rr)/r +
                                    2.0 / (aa + ab) -
                                    3.0 * r / 8.0 / (aa2 + ab2);
                                double ya = -(0.75 + (aa2 + ab2) * 0.25/rr)/r +
                                    2.0 / (aa + ab) -
                                    9.0 * r /16.0/(aa2 + ab2);
                                double ex = xx / r;
                                double ey = yy / r;
                                double ez = zz / r;
                                real0 += (xa - ya) * ex * ex + ya;
                                real1 += (xa - ya) * ex * ey;
                                real2 += (xa - ya) * ex * ez;
                                real3 += (xa - ya) * ey * ey + ya;
                                real4 += (xa - ya) * ey * ez;
                                real5 += (xa - ya) * ez * ez + ya;
                            }
                        }
                    }
                    #pragma vector aligned
                    #pragma simd
                    for (int m = 0; m < nr0; m++) {
                        int mmm = lm0[m];
                        double xx = xx0 + rlx_tbl[mmm];
                        double yy = yy0 + rly_tbl[mmm];
                        double zz = zz0 + rlz_tbl[mmm];
                        double rr = xx * xx + yy * yy + zz * zz;
                        double r = sqrt(rr);
                        double ex = xx / r;
                        double ey = yy / r;
                        double ez = zz / r;
                        double xa;
                        double ya;
                        scalars_ewald_F(xi, r, aa, ab, &xa, &ya);
                        real0 += (xa - ya) * ex * ex + ya;
                        real1 += (xa - ya) * ex * ey;
                        real2 += (xa - ya) * ex * ez;
                        real3 += (xa - ya) * ey * ey + ya;
                        real4 += (xa - ya) * ey * ez;
                        real5 += (xa - ya) * ez * ez + ya;
                    }
                    if (v_real != NULL) {
                        v_real[i0 + 0] += alpha * real0 * f[j0 + 0] +
                                          alpha * real1 * f[j0 + 1] +
                                          alpha * real2 * f[j0 + 2];
                        v_real[i0 + 1] += alpha * real1 * f[j0 + 0] +
                                          alpha * real3 * f[j0 + 1] +
                                          alpha * real4 * f[j0 + 2];
                        v_real[i0 + 2] += alpha * real2 * f[j0 + 0] +
                                          alpha * real4 * f[j0 + 1] +
                                          alpha * real5 * f[j0 + 2];    
                    }
                } // end of if (mat != NULL && mat_real != NULL)
                    
                // recip part
                double recip0 = 0.0;
                double recip1 = 0.0;
                double recip2 = 0.0;
                double recip3 = 0.0;
                double recip4 = 0.0;
                double recip5 = 0.0;
                if (v != NULL || v_recip != NULL) {
                    #pragma vector aligned
                    #pragma simd
                    for (int m = 0; m < nk; m++) {
                        double ex   = ex_tbl[m];
                        double ey   = ey_tbl[m];
                        double ez   = ez_tbl[m];
                        double k1   = k1_tbl[m];
                        double k2   = k2_tbl[m];
                        double k3   = k3_tbl[m];
                        double k    = k_tbl[m];
                        double kexp = kexp_tbl[m];
                        double kk = k * k;
                        double ya = 6.0 * (1.0 - kk * (aa2 + ab2) / 6.0) * kexp;
                        double cf = cos(k1 * xx0 + k2 * yy0 + k3 * zz0);
                        recip0 +=  cf * ya * (1.0 - ex * ex);
                        recip1 += -cf * ya * ex * ey;
                        recip2 += -cf * ya * ex * ez;
                        recip3 +=  cf * ya * (1.0 - ey * ey);
                        recip4 += -cf * ya * ey * ez;
                        recip5 +=  cf * ya * (1.0 - ez * ez);
                    }
                    if (v_recip != NULL) {
                        v_recip[i0 + 0] += alpha * recip0 * f[j0 + 0] +
                                           alpha * recip1 * f[j0 + 1] +
                                           alpha * recip2 * f[j0 + 2];
                        v_recip[i0 + 1] += alpha * recip1 * f[j0 + 0] +
                                           alpha * recip3 * f[j0 + 1] +
                                           alpha * recip4 * f[j0 + 2];
                        v_recip[i0 + 2] += alpha * recip2 * f[j0 + 0] +
                                           alpha * recip4 * f[j0 + 1] +
                                           alpha * recip5 * f[j0 + 2];
                    }
                } // end of if (mat != NULL && mat_recip != NULL)  

                if (v != NULL) {
                    double ewald0 = real0 + recip0;
                    double ewald1 = real1 + recip1;
                    double ewald2 = real2 + recip2;
                    double ewald3 = real3 + recip3;
                    double ewald4 = real4 + recip4;
                    double ewald5 = real5 + recip5;
                    v[i0 + 0] += alpha * ewald0 * f[j0 + 0] +
                                 alpha * ewald1 * f[j0 + 1] +
                                 alpha * ewald2 * f[j0 + 2];
                    v[i0 + 1] += alpha * ewald1 * f[j0 + 0] +
                                 alpha * ewald3 * f[j0 + 1] +
                                 alpha * ewald4 * f[j0 + 2];
                    v[i0 + 2] += alpha * ewald2 * f[j0 + 0] +
                                 alpha * ewald4 * f[j0 + 1] +
                                 alpha * ewald5 * f[j0 + 2];
                } // if (mat != NULL)
            } // for (int j = 0; j < npos; j++) 
        } // for (int i = 0; i < npos; i++)
    } // end of #pragma omp parallel    
}


void NonEwaldKernel(const int npos,
                    const double *pos,
                    const double *rdi,
                    const int ldm,
                    double *mat)
{    
    // compute tile size
    int align_dp = detail::kAlignLen/sizeof(double);
    int tile_size;
    int num_tiles;   
    for (int i = 1; i <= align_dp; i++) {
        if (0 == (3 * i) % align_dp) {
            tile_size = i;
        }
    }
    num_tiles = (npos + tile_size - 1)/tile_size;    
    int64_t num_tiles2 = (int64_t)num_tiles * num_tiles;
    
    #pragma omp parallel                         
    {
        #pragma omp for schedule(dynamic)
        for (int64_t tile = 0; tile < num_tiles2; tile++) {
            int starti = (int)(tile/num_tiles) * tile_size;
            int endi = starti + tile_size;
            endi = endi < npos ? endi : npos;
            int startj = (int)(tile%num_tiles) * tile_size;
            int endj = startj + tile_size;
            endj = endj < npos ? endj : npos;
            if (starti > startj) {
                continue;
            } 
            for (int i = starti; i < endi; i++) {
                const double *pos1 = &(pos[3 * i]);
                double aa = rdi[i];
                double aa2 = aa * aa; 
                int i0 = i * 3;
                for (int j = startj; j < endj; j++) {
                    if (i > j) {
                        continue;
                    }
                    const double *pos2 = &(pos[3 * j]);
                    double ab = rdi[j];
                    double ab2 = ab * ab;
                    int j0 = j * 3;
                    double xa;
                    double ya;
                    double ex = 0.0;
                    double ey = 0.0;
                    double ez = 0.0;
                    if (i != j) {               
                        double xx = pos2[0] - pos1[0];
                        double yy = pos2[1] - pos1[1];
                        double zz = pos2[2] - pos1[2];
                        double rr = xx * xx + yy * yy + zz * zz;
                        double r = sqrt(rr);
                        if (r != 0.0) {
                            ex = xx / r;
                            ey = yy / r;
                            ez = zz / r;
                        }
                        if (r < (aa + ab)) {
                            xa = 2.0 / (aa + ab) - 3.0 * r / 8.0 / (aa2 + ab2);
                            ya = 2.0 / (aa + ab) - 9.0 * r / 16.0 / (aa2 + ab2);
                        } else {
                            xa = (1.5 - (aa2 + ab2) * 0.5 / rr) / r;
                            ya = (0.75 + (aa2 + ab2) * 0.25 / rr) / r;  
                        }
                    } else { // self part
                        xa = ya = 1.0/rdi[i];
                    }
                    double ctmp0 = (xa - ya) * ex * ex + ya;
                    double ctmp1 = (xa - ya) * ex * ey;
                    double ctmp2 = (xa - ya) * ex * ez;
                    double ctmp3 = (xa - ya) * ey * ey + ya;
                    double ctmp4 = (xa - ya) * ey * ez;
                    double ctmp5 = (xa - ya) * ez * ez + ya;
                    mat[(i0 + 0) * ldm + j0 + 0] = ctmp0;
                    mat[(i0 + 0) * ldm + j0 + 1] = ctmp1;
                    mat[(i0 + 0) * ldm + j0 + 2] = ctmp2;
                    mat[(i0 + 1) * ldm + j0 + 0] = ctmp1;
                    mat[(i0 + 1) * ldm + j0 + 1] = ctmp3;
                    mat[(i0 + 1) * ldm + j0 + 2] = ctmp4;
                    mat[(i0 + 2) * ldm + j0 + 0] = ctmp2;
                    mat[(i0 + 2) * ldm + j0 + 1] = ctmp4;
                    mat[(i0 + 2) * ldm + j0 + 2] = ctmp5;
                    if (j0 != i0)
                    {
                        mat[(j0 + 0) * ldm + i0 + 0] = ctmp0;
                        mat[(j0 + 0) * ldm + i0 + 1] = ctmp1;
                        mat[(j0 + 0) * ldm + i0 + 2] = ctmp2;
                        mat[(j0 + 1) * ldm + i0 + 0] = ctmp1;
                        mat[(j0 + 1) * ldm + i0 + 1] = ctmp3;
                        mat[(j0 + 1) * ldm + i0 + 2] = ctmp4;
                        mat[(j0 + 2) * ldm + i0 + 0] = ctmp2;
                        mat[(j0 + 2) * ldm + i0 + 1] = ctmp4;
                        mat[(j0 + 2) * ldm + i0 + 2] = ctmp5;
                    }
                }
            }
        }
    }
}


bool InitRealTable(const double xi, const double box_size,
                   const double tol, EwaldTable *tbl)
{
    tbl->rmax = XitoRmax(xi, tol, box_size, 1.15);
    std::vector<double> rlx_tbl;
    std::vector<double> rly_tbl;
    std::vector<double> rlz_tbl;
    int rmaxi = (int)(tbl->rmax/box_size) + 1;
    double rl = sqrt(3 * box_size * box_size);
    int nr = 0;
    for (int m1 = -rmaxi; m1 <= rmaxi; m1++){
        double rlx = box_size * (double)m1;
        double xx = rlx;
        if (m1 > 0){
            xx -= 0.5 * box_size;
        }
        if (m1 < 0){
            xx += 0.5 * box_size;
        }
        for (int m2 = -rmaxi; m2 <= rmaxi; m2++){
            double rly = box_size * (double) m2;
            double yy = rly;
            if (m2 > 0){
                yy -= 0.5 * box_size;
            }
            if (m2 < 0){
                yy += 0.5 * box_size;
            }
            for (int m3 = -rmaxi; m3 <= rmaxi; m3++){
                double rlz = box_size * (double) m3;
                double zz = rlz;
                if (m3 > 0){
                    zz -= 0.5 * box_size;
                }
                if (m3 < 0){
                    zz += 0.5 * box_size;
                }
                double rr = sqrt(xx * xx + yy * yy + zz * zz);
                if (rr <= tbl->rmax + rl){          
                    rlx_tbl.push_back(rlx);
                    rly_tbl.push_back(rly);
                    rlz_tbl.push_back(rlz);
                    nr++;
                }
            }
        }
    }
    tbl->rlx_tbl = 
        (double *)detail::AlignMalloc(sizeof(double) * nr + detail::kSimdWidth);
    tbl->rly_tbl =
        (double *)detail::AlignMalloc(sizeof(double) * nr + detail::kSimdWidth);
    tbl->rlz_tbl =
        (double *)detail::AlignMalloc(sizeof(double) * nr + detail::kSimdWidth);
    if (NULL == tbl->rlx_tbl ||
        NULL == tbl->rly_tbl ||
        NULL == tbl->rlz_tbl){
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * nr + detail::kSimdWidth);
        return false;
    }
    memcpy(tbl->rlx_tbl, &rlx_tbl[0], sizeof(double) * nr);
    memcpy(tbl->rly_tbl, &rly_tbl[0], sizeof(double) * nr);
    memcpy(tbl->rlz_tbl, &rlz_tbl[0], sizeof(double) * nr);
    tbl->nr = nr;    
    
    return true;
}


bool InitRecipTable(const double xi, const double box_size,
                    const double tol, EwaldTable *tbl)
{
    tbl->kmax = XitoKmax(xi, tol, box_size, 1.15);
    int kmaxi = (int)(tbl->kmax*box_size / (2.0 * M_PI));
    double pivol = M_PI / box_size / box_size / box_size;
    double xi2 = xi * xi;   
    std::vector<double> ex_tbl;
    std::vector<double> ey_tbl;
    std::vector<double> ez_tbl;
    std::vector<double> k1_tbl;
    std::vector<double> k2_tbl;
    std::vector<double> k3_tbl;
    std::vector<double> k_tbl;
    std::vector<double> kexp_tbl;  
    int nk = 0;
    for (double m1 = -kmaxi; m1 <= kmaxi; m1++){
        double k1 = 2.0 * M_PI * (double)m1 / box_size;
        for (double m2 = -kmaxi; m2 <= kmaxi; m2++){
            double k2 = 2.0 * M_PI * (double)m2 / box_size;
            for (double m3 = -kmaxi; m3 <= kmaxi; m3++){
                double k3 = 2.0 * M_PI * (double)m3 / box_size;
                if (m1 != 0.0 || m2 != 0.0 || m3 != 0.0){
                    double kk = k1 * k1 + k2 * k2 + k3 * k3;
                    double k = sqrt(kk);
                    if (k <= tbl->kmax){
                        double k4z = kk / (4.0 * xi2);
                        double k4zplus = 1.0 + k4z * (1.0 + 2.0 * k4z);
                        double kexp = pivol * k4zplus  / kk * exp(-k4z);
                        double ex = k1 / k;
                        double ey = k2 / k;
                        double ez = k3 / k;
                        ex_tbl.push_back(ex);
                        ey_tbl.push_back(ey);
                        ez_tbl.push_back(ez);                       
                        k1_tbl.push_back(k1);
                        k2_tbl.push_back(k2);
                        k3_tbl.push_back(k3);
                        k_tbl.push_back(k);
                        kexp_tbl.push_back(kexp);
                        nk++;
                    }
                }
            }
        }
    }   
    tbl->ex_tbl =
        (double *)detail::AlignMalloc(sizeof(double) * nk + detail::kSimdWidth);
    tbl->ey_tbl =
        (double *)detail::AlignMalloc(sizeof(double) * nk + detail::kSimdWidth);
    tbl->ez_tbl =
        (double *)detail::AlignMalloc(sizeof(double) * nk + detail::kSimdWidth);
    tbl->k1_tbl =
        (double *)detail::AlignMalloc(sizeof(double) * nk + detail::kSimdWidth);
    tbl->k2_tbl =
        (double *)detail::AlignMalloc(sizeof(double) * nk + detail::kSimdWidth);
    tbl->k3_tbl =
        (double *)detail::AlignMalloc(sizeof(double) * nk + detail::kSimdWidth);
    tbl->k_tbl  =
        (double *)detail::AlignMalloc(sizeof(double) * nk + detail::kSimdWidth);
    tbl->kexp_tbl =
        (double *)detail::AlignMalloc(sizeof(double) * nk + detail::kSimdWidth);
    if (NULL == tbl->ex_tbl ||
        NULL == tbl->ey_tbl ||
        NULL == tbl->ez_tbl ||
        NULL == tbl->k1_tbl ||
        NULL == tbl->k2_tbl ||
        NULL == tbl->k3_tbl ||
        NULL == tbl->k_tbl  ||
        NULL == tbl->kexp_tbl){
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * nk + detail::kSimdWidth);
        return false;
    }
    memcpy(tbl->ex_tbl,   &ex_tbl[0],   sizeof(double) * nk);
    memcpy(tbl->ey_tbl,   &ey_tbl[0],   sizeof(double) * nk);
    memcpy(tbl->ez_tbl,   &ez_tbl[0],   sizeof(double) * nk);
    memcpy(tbl->k1_tbl,   &k1_tbl[0],   sizeof(double) * nk);
    memcpy(tbl->k2_tbl,   &k2_tbl[0],   sizeof(double) * nk);
    memcpy(tbl->k3_tbl,   &k3_tbl[0],   sizeof(double) * nk);
    memcpy(tbl->k_tbl,    &k_tbl[0],    sizeof(double) * nk);
    memcpy(tbl->kexp_tbl, &kexp_tbl[0], sizeof(double) * nk);
    tbl->nk = nk;

    return true;
}


bool CreateEwaldTable(const double xi, const double box_size,
                      const double tol, EwaldTable **p_ewald_tbl)
{
    EwaldTable *tbl = new EwaldTable();
    memset(tbl, 0, sizeof(EwaldTable));
    
    InitRealTable(xi, box_size, tol, tbl);
    InitRecipTable(xi, box_size, tol, tbl);
    *p_ewald_tbl = tbl;
    printf("xi = %g nr %d nk %d\n", xi, tbl->nr, tbl->nk);

    return true;
}


void DestroyEwaldTable(EwaldTable *ewald_tbl)
{
    if (ewald_tbl != NULL) {
        detail::AlignFree(ewald_tbl->rlx_tbl);
        detail::AlignFree(ewald_tbl->rly_tbl);
        detail::AlignFree(ewald_tbl->rlz_tbl);
        detail::AlignFree(ewald_tbl->ex_tbl);
        detail::AlignFree(ewald_tbl->ey_tbl);
        detail::AlignFree(ewald_tbl->ez_tbl);
        detail::AlignFree(ewald_tbl->k_tbl);
        detail::AlignFree(ewald_tbl->k1_tbl);
        detail::AlignFree(ewald_tbl->k2_tbl);
        detail::AlignFree(ewald_tbl->k3_tbl);
        detail::AlignFree(ewald_tbl->kexp_tbl);
    }
    delete ewald_tbl;
}


} // namespace stokesdt

} // namespace detail
