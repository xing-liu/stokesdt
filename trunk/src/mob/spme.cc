/**
 * @file   spme.cc
 * @brief  SPME implementation.
 */

#include <mkl.h>
#include <offload.h>
#include <stdint.h>
#include <string.h>

#include "spme.h"


namespace stokesdt {

namespace detail {

#if defined (__MIC__)
const int kSplineLen = 4;
#elif defined (__AVX__) 
const int kSplineLen = 1;
#elif defined (__SSE__)
const int kSplineLen = 1;
#else
const int kSplineLen = 1;
#endif

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
SpmeEngine *spme_mic;
#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)
#define ONCE  alloc_if(1) free_if(1)
#define FREE  alloc_if(0) free_if(1)
#endif


inline double BSpline(int porder, double x)
{
    double spline[porder];
#if 0
    if (1 == n) {
        if (-1 <= x && x < 0) {
            y = 1;
        } else {
            y = 0;
        }
    } else {
        y = (x+1.0)/(n-1.0) * BSpline(n - 1, x) + 
            (n-x-1.0)/(n-1.0) * BSpline(n - 1, x - 1);
    }
#else
    memset(spline, 0, sizeof(double) * porder);
    int p = floor(x) + 1;
    if (p >=0 && p < porder) {
        spline[p] = 1.0;
    }
    for (int i = 2; i <= porder; i++) {
        for (int j = 0; j < porder + 1 - i; j++) {
            spline[j] = (x - j + 1.0) * spline[j] +
                (i - x + j - 1.0) * spline[j + 1];
            spline[j] /= i - 1.0;
        }
    }
#endif
    
    return spline[0];
}


static bool ComputeSpline(const double xi,
                          const int dim,
                          const int porder,
                          const double box_size,
                          const int ld1,
                          const int ld2,
                          double **p_map,
                          double **p_lm2)
{
    double splineval[porder - 1];
    for (int i = 0; i < porder - 1; i++) {
        splineval[i] = BSpline(porder, i);
    }
    
    double *map = (double *)AlignMalloc(sizeof(double) * dim + kSimdWidth);
    if (map == NULL) {
        return false;
    }
    double *lm2 = (double *)AlignMalloc(sizeof(double) * ld2/2 * dim);
    if (lm2 == NULL) {
        return false;    
    }
    MKL_Complex16 *b = 
        (MKL_Complex16 *)AlignMalloc(sizeof(MKL_Complex16)*dim + kSimdWidth);
    if (b == NULL) {
        return false;    
    }

    // compute map
    int x = 0;
    for (int i = 0; i <= dim/2; i++) {
        map[x] = i * 2.0 * M_PI / box_size;
        x++;
    }
    for (int i = -(dim/2-1); i <= -1; i++) {
        map[x] = i * 2.0 * M_PI / box_size;
        x++;
    }

    // compute b
    for (int x = 0; x < dim; x++) {
        double kx = map[x];
        MKL_Complex16 den;
        MKL_Complex16 btmp;
        den.real = 0.0;
        den.imag = 0.0;
        for (int i = 0; i < porder-1; i++) {
            double phi = i * kx * box_size / dim;
            den.real += splineval[i] * cos(phi);
            den.imag += splineval[i] * sin(phi);
        }
        double phi0 = box_size * kx * (porder-1) / dim;
        btmp.real = cos(phi0);
        btmp.imag = sin(phi0);        
        double m2 = den.real*den.real + den.imag*den.imag;
        b[x].real = (btmp.real * den.real + btmp.imag * den.imag)/m2;
        b[x].imag = (btmp.imag * den.real - btmp.real * den.imag)/m2;
    }

    // compute lm2
    for (int j = 0; j < dim * dim; j++) {
        int z = j/dim;
        int y = j%dim;
        double kz = map[z];
        double ky = map[y];      
        for (int x = 0; x <= dim/2; x+=kSplineLen) {                   
            for (int i = 0; i < kSplineLen; i++) {
                MKL_Complex16 btmp;
                MKL_Complex16 bbb;
                double kx = map[x + i];
                double k = sqrt(kx*kx + ky*ky + kz*kz);
                double m2 = ScalarRecip(k, xi);
                btmp.real =
                    b[x + i].real * b[y].real - b[x + i].imag * b[y].imag;
                btmp.imag =
                    b[x + i].real * b[y].imag + b[x + i].imag * b[y].real;
                bbb.real = btmp.real * b[z].real - btmp.imag * b[z].imag;
                bbb.imag = btmp.real * b[z].imag + btmp.imag * b[z].real;
                m2 = m2 * (bbb.real*bbb.real + bbb.imag*bbb.imag);
                lm2[j * ld1/2 + x + i] = m2;
            }
        }
    }
    
    *p_map = map;
    *p_lm2 = lm2;
    AlignFree(b);

    return true;
}


static void Interpolate(const int npos,
                        const double *rdi,
                        const int porder3,
                        const double *P,
                        const int ldP,
                        const int *ind,
                        const int ldind,                  
                        const double alpha,
                        const double *grid,
                        const int ld3,
                        const double beta,
                        double *vels )
{
    const int align_dp = kAlignLen/sizeof(double);

    if (rdi != NULL) {
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < npos; i+=align_dp) {
            int endi = i + align_dp > npos ? npos : i + align_dp;
            for (int k = i; k < endi; k++) {
                double alpha0 = rdi[k] * alpha;
                InterpolateKernel(porder3, &(P[k*ldP]), &(ind[k*ldind]),
                                  alpha0, grid, ld3, beta, &(vels[3*k]));
            }
        }
    } else {
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < npos; i+=align_dp) {
            int endi = i + align_dp > npos ? npos : i + align_dp;
            for (int k = i; k < endi; k++) {
                InterpolateKernel(porder3, &(P[k*ldP]), &(ind[k*ldind]),
                                  alpha, grid, ld3, beta, &(vels[3*k]));
            }
        }    
    } // if (rdi != NULL)
}


static void ApplyInfluence(const int flag,
                           const int dim,
                           const double *map,
                           const double *lm2,                     
                           const int ld1,
                           const int ld2,
                           const int ld3,
                           double *grid)
{    
    int ld1c = ld1/2;
    int ld2c = ld2/2;
    
    // grid(x,y,z,:) = B(:,:,x,y,z)*squeeze(grid(x,y,z,:));    
    #pragma omp parallel default(none)\
                         shared (ld1c, ld2c, ld3, dim, grid, map, lm2, flag)
    {
        __declspec(align(kAlignLen)) double B0[dim/2 + 1 + kSplineLen];
        __declspec(align(kAlignLen)) double B1[dim/2 + 1 + kSplineLen];
        __declspec(align(kAlignLen)) double B2[dim/2 + 1 + kSplineLen];
        __declspec(align(kAlignLen)) double B3[dim/2 + 1 + kSplineLen];
        __declspec(align(kAlignLen)) double B4[dim/2 + 1 + kSplineLen];
        __declspec(align(kAlignLen)) double B5[dim/2 + 1 + kSplineLen];

        #pragma omp for schedule(static)
        for (int j = 0; j < dim * dim; j++) {
            int z = j/dim;
            int y = j%dim;
            double kz = map[z];
            double ky = map[y];
            const double *mm = &(lm2[j * ld1c]);
            double kyz = ky*ky + kz*kz;
            if (flag == 0) {
                #pragma simd
                for (int x = 0; x <= dim/2; x++) {
                    double kx = map[x];
                    double m2 = mm[x];
                    double kk = kx*kx + kyz;
                    double kk_1 = 1.0/kk; 
                    m2 =  m2 * (1.0 - kk/3.0);
                    double exx = kx * kx * kk_1;
                    double eyy = ky * ky * kk_1;
                    double ezz = kz * kz * kk_1;
                    double exy = kx * ky * kk_1;
                    double eyz = ky * kz * kk_1;
                    double exz = kx * kz * kk_1;
                    B0[x] = m2 - m2 * exx; // 0
                    B1[x] = -m2 * exy;     // 1 and 3
                    B2[x] = -m2 * exz;     // 2 and 6
                    B3[x] = m2 - m2 * eyy; // 4
                    B4[x] = -m2 * eyz;     // 5 and 7
                    B5[x] = m2 - m2 * ezz; // 8
                }
            } else if (flag == 1) {
                #pragma simd
                for (int x = 0; x <= dim/2; x++) {
                    double kx = map[x];
                    double m2 = mm[x];
                    double kk = kx*kx + kyz;
                    double kk_1 = 1.0/kk;
                    double exx = kx * kx * kk_1;
                    double eyy = ky * ky * kk_1;
                    double ezz = kz * kz * kk_1;
                    double exy = kx * ky * kk_1;
                    double eyz = ky * kz * kk_1;
                    double exz = kx * kz * kk_1;
                    B0[x] = m2 - m2 * exx; // 0
                    B1[x] = -m2 * exy;     // 1 and 3
                    B2[x] = -m2 * exz;     // 2 and 6
                    B3[x] = m2 - m2 * eyy; // 4
                    B4[x] = -m2 * eyz;     // 5 and 7
                    B5[x] = m2 - m2 * ezz; // 8
                }            
            } else if (flag == 2) {
                #pragma simd
                for (int x = 0; x <= dim/2; x++) {
                    double kx = map[x];
                    double m2 = mm[x];
                    double kk = kx*kx + kyz;
                    double kk_1 = 1.0/kk;
                    m2 =  m2 * (-kk/6.0);                    
                    double exx = kx * kx * kk_1;
                    double eyy = ky * ky * kk_1;
                    double ezz = kz * kz * kk_1;
                    double exy = kx * ky * kk_1;
                    double eyz = ky * kz * kk_1;
                    double exz = kx * kz * kk_1;
                    B0[x] = m2 - m2 * exx; // 0
                    B1[x] = -m2 * exy;     // 1 and 3
                    B2[x] = -m2 * exz;     // 2 and 6
                    B3[x] = m2 - m2 * eyy; // 4
                    B4[x] = -m2 * eyz;     // 5 and 7
                    B5[x] = m2 - m2 * ezz; // 8
                }             
            }

            for (int x = 0; x <= dim/2; x+=kSplineLen) {
                int i = z * ld2c + y * ld1c + x;
                InfluenceKernel(&(B0[x]), &(B1[x]), &(B2[x]),
                                &(B3[x]), &(B4[x]), &(B5[x]),
                                ld3, &(grid[2*i]));
            } 
        }
    }
    grid[0 * ld3 + 0] = 0.0;
    grid[0 * ld3 + 1] = 0.0;
    grid[1 * ld3 + 0] = 0.0;
    grid[1 * ld3 + 1] = 0.0;
    grid[2 * ld3 + 0] = 0.0;
    grid[2 * ld3 + 1] = 0.0;
}


static void Spread(const double *rdi,
                   const int nb,
                   const int *head,
                   const int *next,
                   const int porder3,
                   const double *P,
                   const int ldP,
                   const int *ind,
                   const int ldind,
                   const double *forces,
                   const int ld3,
                   double *grid)
{    
    int nb3 = nb * nb * nb;
    int n8 = nb3 / 8;
   
    #pragma omp parallel for
    for (int i = 0; i < ld3 * 3; i++) {
        grid[i] = 0.0;
    }

    if (rdi != NULL) {
        for (int set = 0; set < 8; set++) {
            #pragma omp parallel for
            for (int k = n8 * set; k < (set + 1) * n8; k++) {
                int i = head[k];
                while (i != -1) {
                    double force0[3];
                    force0[0] = forces[i * 3 + 0] * rdi[i];
                    force0[1] = forces[i * 3 + 1] * rdi[i];
                    force0[2] = forces[i * 3 + 2] * rdi[i];
                    SpreadKernel(porder3, &(P[i * ldP]), &(ind[i * ldind]),
                                 force0, ld3, grid);
                    i = next[i];
                }
            }
        }
    } else {
        for (int set = 0; set < 8; set++) {
            #pragma omp parallel for
            for (int k = n8 * set; k < (set + 1) * n8; k++) {
                int i = head[k];
                while (i != -1) {
                    SpreadKernel(porder3, &(P[i * ldP]), &(ind[i * ldind]),
                                 &(forces[i * 3]), ld3, grid);
                    i = next[i];
                }
            }
        }    
    } // if (rdi != NULL)
}


static void UpdateSpmeEngine_(const double *pos,
                              SpmeEngine *spme)
{
    int npos = spme->npos;
    double pidx[npos];
    #pragma omp parallel default(none)\
                         shared (npos, pos, pidx, spme)
    {
        int porder = spme->porder;
        int dim = spme->dim;
//        const double *W = tab_W[porder - 1];
        int porder2 = porder * porder;
        int ldP = spme->ldP;
        double *P = spme->P;
        int ldind = spme->ldind;
        int *ind = spme->ind;
        int sizeb = spme->sizeb;
        int ld1 = spme->ld1;
        int ld2 = spme->ld2;
        int nb = spme->nb;
        int nb2 = nb * nb;
        int nb3 = nb2 * nb;
        double box_size = spme->box_size;
        int *bidx = spme->bidx;
        int *head = spme->head;
        int *next = spme->next;
    #if 0        
        __declspec(align(kAlignLen)) double p1[porder];
        __declspec(align(kAlignLen)) double p2[porder];
        __declspec(align(kAlignLen)) double p3[porder];
    #endif
        __declspec(align(kAlignLen)) double q1[porder];
        __declspec(align(kAlignLen)) double q2[porder];
        __declspec(align(kAlignLen)) double q3[porder];        
        double g[3];
        double ground[3];
        #pragma omp for schedule(static)
        for (int i = 0; i < npos; i++) {
            double cx = fmod(pos[3 * i + 0], box_size);
            double cy = fmod(pos[3 * i + 1], box_size);
            double cz = fmod(pos[3 * i + 2], box_size);
            cx = (cx >= 0.0 ? cx : box_size + cx);
            cx = (cx >= 0.0 ? cx : box_size + cx);
            cx = (cx >= 0.0 ? cx : box_size + cx);            
            cx = pos[i * 3 + 0]/box_size * dim;
            cy = pos[i * 3 + 1]/box_size * dim;
            cz = pos[i * 3 + 2]/box_size * dim;
            // independent set
            int bx = (int)(cx/sizeb);
            int by = (int)(cy/sizeb);
            int bz = (int)(cz/sizeb);
            pidx[i] = bidx[bz * nb2 + by * nb + bx];         
            // FFT mesh
            g[0] = floor(cx);
            g[1] = floor(cy);
            g[2] = floor(cz);
            int indx;
            int indy;
            int indz;
            if (porder%2 == 0) {
                indx = (int)(g[0] - porder/2 + 1);
                indx = (indx + dim) % dim;
                indy = (int)(g[1] - porder/2 + 1);
                indy = (indy + dim) % dim;
                indz = (int)(g[2] - porder/2 + 1);
                indz = (indz + dim) % dim;
            } else {
                ground[0] = floor(cx + 0.5);
                ground[1] = floor(cy + 0.5);
                ground[2] = floor(cz + 0.5);
                int temp = (porder - 1)/2;
                indx = (int)(g[0] - temp);
                indx = (indx + dim) % dim;
                indy = (int)(g[1] - temp);
                indy = (indy + dim) % dim;
                indz = (int)(g[2] - temp);
                indz = (indz + dim) % dim;
            }
        #if 0
            g[0] = cx - g[0] - 0.5;
            g[1] = cy - g[1] - 0.5;
            g[2] = cz - g[2] - 0.5;
        
            p1[0] = 1.0;
            p2[0] = 1.0;
            p3[0] = 1.0;
            for (int x = 1; x < porder; x++) {
                p1[x] = g[0] * p1[x - 1];
                p2[x] = g[1] * p2[x - 1];
                p3[x] = g[2] * p3[x - 1];
            }
        
            for (int x = 0; x < porder; x++) {
                q1[x] = 0.0;
                q2[x] = 0.0;
                q3[x] = 0.0;           
                for (int y = 0; y < porder; y++) {
                    q1[x] += W[x * porder + y] * p1[y];
                    q2[x] += W[x * porder + y] * p2[y];
                    q3[x] += W[x * porder + y] * p3[y];
                }
            }
        #else
            for (int x = 0; x < porder; x++) {
                q1[x] = BSpline(porder, g[0] - cx + x);
                q2[x] = BSpline(porder, g[1] - cy + x);
                q3[x] = BSpline(porder, g[2] - cz + x);
            }
        #endif
            for (int z = 0; z < porder; z++) {
                int zz = (z + indz) % dim;
                for (int y = 0; y < porder; y++) {
                    int yy = (y + indy) % dim;
                    for (int x = 0; x < porder; x++) {
                        int xx = (x + indx) % dim;
                        P[i * ldP + z * porder2 + y * porder + x]
                            = q1[x] * q2[y] * q3[z];
                        ind[i * ldind + z * porder2 + y * porder + x]
                            = zz * ld2 + yy * ld1 + xx;
                    }
                }
            }
        }

        int tid = omp_get_thread_num ();
        int nthreads = omp_get_num_threads();
        int startidx = (nb3 + nthreads - 1)/nthreads * tid;
        int endidx = (nb3 + nthreads - 1)/nthreads * (tid + 1);        
        // init myhead
        #pragma omp for
        for (int i = 0; i < nb3; i++) {
            head[i] = -1;
        }
        // compute myhead 
        for (int i = 0; i < npos; i++) {
            int idx = pidx[i];
            if (idx >= startidx && idx < endidx) {
                int h = head[idx];
                head[idx] = i;
                next[i] = h;
            }
        }
    } /* #pragma omp parallel */
}


static bool CreateSpmeEngine_(const int npos,
                              const double *rdi,
                              const double box_size,
                              const double xi,
                              const int dim,
                              const int porder,
                              SpmeEngine **p_spme)
{    
    SpmeEngine *spme = new SpmeEngine;
    spme->dim = dim;
    spme->porder = porder;    
    spme->npos = npos;
    spme->box_size = box_size;
    spme->xi = xi;
    if (rdi != NULL) {
        spme->rdi2 = (double *)AlignMalloc(sizeof(double) * npos);
        if (NULL == spme->rdi2) {
            return false;
        }
        for (int i = 0; i < npos; i++) {
            spme->rdi2[i] = rdi[i] * rdi[i];
        }
    } else {
        spme->rdi2 = NULL;
    }
        
    // init fft
    int pad_dim  = FFTPadLen(dim, sizeof(double));
    int pad_dim2 = (dim/2 + 1)*2;
    pad_dim2 = FFTPadLen(pad_dim2, sizeof(double));
    spme->ld1 = pad_dim2;
    spme->ld2 = pad_dim * pad_dim2;
    spme->ld3 = dim * pad_dim * pad_dim2;  
    spme->grid = (double *)AlignMalloc(sizeof(double) * spme->ld3 * 3);
    if (NULL == spme->grid) {
        return false;
    }
    MKL_LONG size[3];
    MKL_LONG ld_fw[4];
    MKL_LONG ld_bw[4];
    size[0] = dim;
    size[1] = dim;
    size[2] = dim;
    ld_fw[0] = 0;
    ld_fw[1] = spme->ld2;
    ld_fw[2] = spme->ld1;
    ld_fw[3] = 1;    
    ld_bw[0] = 0;
    ld_bw[1] = spme->ld2/2;
    ld_bw[2] = spme->ld1/2;
    ld_bw[3] = 1;

    // r2c FFT
    MKL_LONG ret = DftiCreateDescriptor(&(spme->fwhandle),
                                        DFTI_DOUBLE, DFTI_REAL, 3, size);
    if (ret != 0) {
        return false;            
    }

    ret = DftiSetValue(spme->fwhandle, DFTI_INPUT_STRIDES, ld_fw);
    if (ret != 0) {
        return false;            
    }
    ret = DftiSetValue(spme->fwhandle, DFTI_OUTPUT_STRIDES, ld_bw);
    if (ret != 0) {
        return false;            
    }
    ret = DftiSetValue(spme->fwhandle,
                       DFTI_CONJUGATE_EVEN_STORAGE,
                       DFTI_COMPLEX_COMPLEX);
    if (ret != 0) {
        return false;            
    }

    ret = DftiCommitDescriptor(spme->fwhandle);
    if (ret != 0) {
        return false;            
    }
    
    // c2r FFT
    ret = DftiCreateDescriptor(&(spme->bwhandle),
                               DFTI_DOUBLE, DFTI_REAL, 3, size);
    if (ret != 0) {
        return false;            
    }
    ret = DftiSetValue(spme->bwhandle, DFTI_INPUT_STRIDES, ld_bw); 
    if (ret != 0) {
        return false;            
    }  
    ret = DftiSetValue(spme->bwhandle, DFTI_OUTPUT_STRIDES,ld_fw);
    if (ret != 0) {
        return false;            
    }
    ret = DftiSetValue(spme->bwhandle,
                       DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (ret != 0) {
        return false;            
    }
    ret = DftiSetValue(spme->bwhandle,
                       DFTI_BACKWARD_SCALE, 1.0/((double)dim*dim*dim));
    if (ret != 0) {
        return false;            
    }
    DftiCommitDescriptor(spme->bwhandle);
    if (ret != 0) {
        return false;            
    }

    // compute spline
    ComputeSpline(xi, dim, porder, box_size, spme->ld1, spme->ld2,
                  &(spme->map), &(spme->lm2));
           
    // sparse P
    int porder2 = porder * porder;
    int porder3 = porder2 * porder;
    spme->ldP = PadLen(porder3, sizeof(double));
    spme->P = (double *)AlignMalloc(sizeof(double) * spme->ldP * npos);   
    spme->ldind = PadLen(porder3, sizeof(int));
    spme->ind = (int *)AlignMalloc(sizeof(int) * spme->ldind * npos);
    if (NULL == spme->P || NULL == spme->ind) {
        return false;
    }
    
    // init independent set
    spme->sizeb = spme->porder;
    while (spme->dim%spme->sizeb != 0)
    {
        spme->sizeb++;        
    }
    int nb = spme->dim/spme->sizeb;
    spme->nb = nb;
    if (nb < 4) {
        return false;
    }
    spme->next = (int *)malloc(sizeof(int) * npos);
    if (NULL == spme->next) {
        return false;
    }
    spme->bidx = (int *)malloc(sizeof(int) * nb * nb * nb);
    if (NULL == spme->bidx) {
        return false;
    }
    spme->head = (int *)malloc(sizeof(int) * nb * nb * nb);
    if (NULL == spme->head) {
        return false;
    }
    for (int i = 0; i < nb * nb * nb; i++) {
        spme->head[i] = -1;
    }

    int iset = 0;
    for (int sx = 0; sx < 2; sx++) {
        for (int sy = 0; sy < 2; sy++) {
            for (int sz = 0; sz < 2; sz++) {
                for (int x = sx; x < nb; x+=2) {
                    for (int y = sy; y < nb; y+=2) {
                        for (int z = sz; z < nb; z+=2) {
                            int idx = z * nb * nb + y * nb + x;
                            spme->bidx[idx] = iset;
                            iset++;
                        }
                    }              
                }
            }
        }
    }

    // NUMA
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < spme->ld3; i++) {
        spme->grid[0 * spme->ld3 + i] = 0.0;
        spme->grid[1 * spme->ld3 + i] = 0.0;
        spme->grid[2 * spme->ld3 + i] = 0.0;           
    }

    *p_spme = spme;

    return true;
}


static void DestroySpmeEngine_(SpmeEngine *spme)
{
    if (spme != NULL) {
        DftiFreeDescriptor(&(spme->fwhandle));
        DftiFreeDescriptor(&(spme->bwhandle));        
        AlignFree(spme->map);
        AlignFree(spme->lm2);
        AlignFree(spme->grid);
        AlignFree(spme->P);
        AlignFree(spme->ind);
        free(spme->head);
        free(spme->bidx);
        free(spme->next);
        if (spme->rdi2 != NULL) {
            AlignFree(spme->rdi2);
        }
    }
    delete spme;
}


static void ComputeSpmeRecip_(const SpmeEngine *spme,
                              const int nrhs,
                              const double alpha,                              
                              const int ldin,
                              const double *vec_in,
                              const double beta,
                              const int ldout,
                              double *vec_out)
{
#ifdef __INTEL_OFFLOAD
    int idev = _Offload_get_device_number();
#endif
    int npos = spme->npos;
    double box_size = spme->box_size;
    int ld1 = spme->ld1;
    int ld2 = spme->ld2;
    int ld3 = spme->ld3;
    int dim = spme->dim;
    int dim2 = dim * dim;
    int dim3 = dim2 * dim;
    int porder = spme->porder;
    int porder3 = porder * porder * porder;
    double *P = spme->P;
    int *ind = spme->ind;
    int ldP = spme->ldP;
    int ldind = spme->ldind;
    double *grid = spme->grid;
    double alpha0 = alpha * (double)dim3/(box_size*box_size*box_size);

    if (spme->rdi2 == NULL) {
        for (int irhs = 0; irhs < nrhs; irhs++) {
            const double *vin = &(vec_in[ldin * irhs]);
            double *vout = &(vec_out[ldout * irhs]);

            // spread
            Spread(NULL, spme->nb, spme->head, spme->next,
                   porder3, P, ldP, ind, ldind,
                   vin, ld3, grid);

            // forward fft
            DftiComputeForward(spme->fwhandle, &(grid[0 * ld3]));
            DftiComputeForward(spme->fwhandle, &(grid[1 * ld3]));
            DftiComputeForward(spme->fwhandle, &(grid[2 * ld3]));

            // apply influence
            ApplyInfluence(0, dim, spme->map, spme->lm2,
                           ld1, ld2, ld3, grid);

            // backward fft
            DftiComputeBackward(spme->bwhandle, &(grid[0 * ld3]));
            DftiComputeBackward(spme->bwhandle, &(grid[1 * ld3]));
            DftiComputeBackward(spme->bwhandle, &(grid[2 * ld3]));

            // interpolate
            Interpolate(npos, NULL, porder3, P, ldP, ind, ldind,
                        alpha0, grid, ld3, beta, vout);
        }
    } else {
        for (int irhs = 0; irhs < nrhs; irhs++) {
            const double *vin = &(vec_in[ldin * irhs]);
            double *vout = &(vec_out[ldout * irhs]);

            // 1st call
            // spread
            Spread(NULL, spme->nb, spme->head, spme->next,
                   porder3, P, ldP, ind, ldind,
                   vin, ld3, grid);
            // forward fft
            DftiComputeForward (spme->fwhandle, &(grid[0 * ld3]));
            DftiComputeForward (spme->fwhandle, &(grid[1 * ld3]));
            DftiComputeForward (spme->fwhandle, &(grid[2 * ld3]));
            // apply influence
            ApplyInfluence(1, dim, spme->map, spme->lm2,
                           ld1, ld2, ld3, grid);
            // backward fft
            DftiComputeBackward (spme->bwhandle, &(grid[0 * ld3]));
            DftiComputeBackward (spme->bwhandle, &(grid[1 * ld3]));
            DftiComputeBackward (spme->bwhandle, &(grid[2 * ld3]));
            // interpolate
            Interpolate(npos, NULL, porder3, P, ldP, ind, ldind,
                        alpha0, grid, ld3, beta, vout);

            // 2nd call
            // spread
            Spread(spme->rdi2, spme->nb, spme->head, spme->next,
                   porder3, P, ldP, ind, ldind,
                   vin, ld3, grid);
            // forward fft
            DftiComputeForward (spme->fwhandle, &(grid[0 * ld3]));
            DftiComputeForward (spme->fwhandle, &(grid[1 * ld3]));
            DftiComputeForward (spme->fwhandle, &(grid[2 * ld3]));
            // apply influence
            ApplyInfluence(2, dim, spme->map, spme->lm2,
                           ld1, ld2, ld3, grid);
            // backward fft
            DftiComputeBackward (spme->bwhandle, &(grid[0 * ld3]));
            DftiComputeBackward (spme->bwhandle, &(grid[1 * ld3]));
            DftiComputeBackward (spme->bwhandle, &(grid[2 * ld3]));
            // interpolate
            Interpolate(npos, NULL, porder3, P, ldP, ind, ldind,
                        alpha0, grid, ld3, 1.0, vout);


            // 3rd call
            // spread
            Spread(NULL, spme->nb, spme->head, spme->next,
                   porder3, P, ldP, ind, ldind,
                   vin, ld3, grid);
            // forward fft
            DftiComputeForward (spme->fwhandle, &(grid[0 * ld3]));
            DftiComputeForward (spme->fwhandle, &(grid[1 * ld3]));
            DftiComputeForward (spme->fwhandle, &(grid[2 * ld3]));
            // apply influence
            ApplyInfluence(2, dim, spme->map, spme->lm2,
                           ld1, ld2, ld3, grid);
            // backward fft
            DftiComputeBackward (spme->bwhandle, &(grid[0 * ld3]));
            DftiComputeBackward (spme->bwhandle, &(grid[1 * ld3]));
            DftiComputeBackward (spme->bwhandle, &(grid[2 * ld3]));
            // interpolate
            Interpolate(npos, spme->rdi2, porder3, P, ldP, ind, ldind,
                        alpha0, grid, ld3, 1.0, vout);
            
        }
    }
}


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


bool CreateSpmeEngine(const int npos,
                      const double *rdi,
                      const double box_size,
                      const double xi,
                      const int dim,
                      const int porder,
                      SpmeEngine **p_spme)
{
    bool ret;
    ret = CreateSpmeEngine_(npos, rdi, box_size, xi, dim, porder, p_spme);
    if (!ret) {
        return false;
    }
    
#ifdef __INTEL_OFFLOAD
    int nbcpu = 1;
    int nbmic = 1;
    int mic_numdevs = _Offload_number_of_devices();;
    int nbs = nbmic * mic_numdevs + nbcpu;
    p_spme[0]->nbs = nbs;
    p_spme[0]->nbcpu = nbcpu;
    p_spme[0]->nbmic = nbmic;
    p_spme[0]->mic_numdevs = mic_numdevs;

    if (NULL == rdi) {
        for (int i = 0; i < mic_numdevs; i++) {
            #pragma offload target(mic: i)\
                    in(xi, dim, porder, npos, box_size)\
                    in(rdi :length(npos) ONCE)\
                    out(ret) nocopy(spme_mic)
            ret = CreateSpmeEngine_(npos, rdi, box_size,
                                    xi, dim, porder, &spme_mic);
            if (!ret) {
                return false;
            }
        }
    } else {
        for (int i = 0; i < mic_numdevs; i++) {
            #pragma offload target(mic: i)\
                    in(xi, dim, porder, npos, box_size)\
                    in(rdi :length(npos) ONCE)\
                    out(ret) nocopy(spme_mic)
            ret = CreateSpmeEngine_(npos, rdi, box_size,
                                    xi, dim, porder, &spme_mic);
            if (!ret) {
                return false;
            }
        }    
    }
#endif

    return true;
}


void DestroySpmeEngine(SpmeEngine *spme)
{
#ifdef __INTEL_OFFLOAD
    for (int i = 0; i < spme->mic_numdevs; i++) {
        #pragma offload target(mic: i)\
                nocopy(spme_mic)
        DestroySpmeEngine_(spme_mic);
    }
#endif
    DestroySpmeEngine_(spme);
}

    
void ComputeSpmeRecip(const SpmeEngine *spme,
                      const int nrhs,
                      const double alpha,                              
                      const int ldin,
                      const double *vec_in,
                      const double beta,
                      const int ldout,
                      double *vec_out)
{  
#ifdef __INTEL_OFFLOAD
    int nbmic = spme->nbmic;
    int nbcpu = spme->nbcpu;
    int nbs = spme->nbs;
    int nm = spme->npos * 3;
    int mic_numdevs = spme->mic_numdevs;
    int niters = (nrhs + nbs - 1) / nbs;
    int idx = 0;
    for (int i = 0; i < niters; i++) {
        int workdevs = 0;
        if (nbmic != 0 && nrhs != 1) {
            workdevs = 0;
            for (int j = 0; j < mic_numdevs; j++) {
                int size = idx + nbmic > nrhs ? nrhs - idx : nbmic;   
                const double *f = &(vec_in[idx * ldin]);
                double *v = &(vec_out[idx * ldout]);

                #pragma offload target(mic:j) signal(j)\
                            in(alpha, beta, ldin, ldout, size)\
                            in(f :length(ldin*size) ONCE)\
                            inout(v :length(ldout*size) ONCE)\
                            nocopy(spme_mic)
                ComputeSpmeRecip_(spme_mic, nrhs, alpha, ldin, f,
                                  beta, ldout, v);                
                workdevs++;
                idx += size;
                if (idx == nrhs) {
                    break;
                }
            }
        }

        // CPU work concurrently
        if (idx != nrhs && nbcpu != 0) {
            int size = idx + nbcpu > nrhs ? nrhs - idx : nbcpu;
            ComputeSpmeRecip_(spme, nrhs, alpha, ldin, &(vec_in[idx * ldin]),
                              beta, ldout, &(vec_out[idx * ldout]));
            idx += size;
        }
         
        // wait for complete
        for (int j = 0; j < workdevs; j++) {
            #pragma offload_wait target(mic:j) wait(j)
        }
        
        if (idx == nrhs) {
            break;
        }
    }
#else
    ComputeSpmeRecip_(spme, nrhs, alpha, ldin,
                      vec_in, beta, ldout, vec_out);
#endif // #ifdef __INTEL_OFFLOAD
}


void UpdateSpmeEngine(const double *pos, SpmeEngine *spme)
{
#ifdef __INTEL_OFFLOAD
    int npos = spme->npos;
    for (int i = 0; i < spme->mic_numdevs; i++) {
        #pragma offload target(mic: i)\
                in(pos :length(3*npos) ONCE)\ 
                nocopy(spme_mic)
        UpdateSpmeEngine_(pos, spme_mic);
    }
#endif
    UpdateSpmeEngine_(pos, spme);
}

} // namespace stokesdt

} // namespace detail