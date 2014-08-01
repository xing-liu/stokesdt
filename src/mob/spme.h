/**
 * @file   spme.h
 * @brief  SPME engine definition
 */

#ifndef SPME_H_
#define SPME_H_


#include <mkl.h>
#include <cmath>
#include <x86intrin.h>
#include "common.h"


namespace stokesdt {

namespace detail {
    

/** \struct  SpmeEngine
 *  @brief   SPME engine
 */
typedef struct SpmeEngine {
#ifdef __INTEL_OFFLOAD
    int mic_numdevs;
    int nbmic;
    int nbcpu;
    int nbs;
#endif
    /// the Ewald parameter
    double xi;
    /// the number of particles
    int npos;
    /// the dimension of the simulation box
    double box_size;

    /// the interpolation order
    int porder; 
    /// the matrix P (in sparse format), \c porder^3 * np
    double *P;
    /// the leading dimension of the matrix P, \c porder^3
    int ldP;
    /// the array of column indices of the matrix P
    int *ind;
    /// the leading dimension of ind
    int ldind;

    /// the pointer to the influence map
    double *map;
    /// the array of Ewald factors
    double *lm2;    
    
    /// the dimension of independent sets    
    int nb;
    /// the size of independent block
    int sizeb;
    /// the cell list of particles, head, \c nb * nb * nb
    int *head;
    /// the cell list of particles, bidx, \c nb * nb * nb
    int *bidx;
    /// the cell list of particles, next, \c np
    int *next;

    /// the forward FFT handle
    DFTI_DESCRIPTOR_HANDLE bwhandle;
    /// the backward FFT handle
    DFTI_DESCRIPTOR_HANDLE fwhandle;    
    /// the FFT mesh
    double *grid;
    /// the dimension of the FFT mesh
    int dim;
    /// the first leading dimension of the FFT mesh
    int ld1;
    /// the second leading dimension of the FFT mesh
    int ld2;
    /// the third leading dimension of the FFT mesh
    int ld3;
} SpmeEngine;


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

const double tab_splines[6][5] = {
  {0.0},
  {0.0},
  {0.0},
  {1.0/6.0, 4.0/6.0, 1.0/6.0},
  {0.0},
  {1.0/120.0, 26.0/120.0, 66.0/120.0, 26.0/120.0, 1.0/120.0}
};


const double tab_W[6][6*6]= {
  {0.0},
  {0.0},
  {0.0},
  { 1.0/48.0,  -6.0/48.0,  12.0/48.0,  -8.0/48.0,
   23.0/48.0, -30.0/48.0, -12.0/48.0,  24.0/48.0,
   23.0/48.0,  30.0/48.0, -12.0/48.0, -24.0/48.0,
    1.0/48.0,   6.0/48.0,  12.0/48.0,   8.0/48.0},
  {0.0},
  {   1.0/3840.0,   -10.0/3840.0,   40.0/3840.0,   -80.0/3840.0,   80.0/3840.0,  -32.0/3840.0,
    237.0/3840.0,  -750.0/3840.0,  840.0/3840.0,  -240.0/3840.0, -240.0/3840.0,  160.0/3840.0,
   1682.0/3840.0, -1540.0/3840.0, -880.0/3840.0,  1120.0/3840.0,  160.0/3840.0, -320.0/3840.0,
   1682.0/3840.0,  1540.0/3840.0, -880.0/3840.0, -1120.0/3840.0,  160.0/3840.0,  320.0/3840.0,
    237.0/3840.0,   750.0/3840.0,  840.0/3840.0,   240.0/3840.0, -240.0/3840.0, -160.0/3840.0,
      1.0/3840.0,    10.0/3840.0,   40.0/3840.0,    80.0/3840.0,   80.0/3840.0,   32.0/3840.0}  
};


inline double ScalarRecip(const double k,
                          const double xi,
                          const double aa,
                          const double ab)
{
    double kx = k / xi;
    double kx2 = kx * kx;
    double k2 = k * k;
    double aa2 = aa*aa;
    double ab2 = ab*ab;
    double v = 6.0 * M_PI * exp(-kx2/4.0)
        * (1.0 - k2*(aa2+ab2)/6.0) / k2 * (1.0 + kx2 * (1.0/4.0 + kx2/8.0));
    return v;
}


inline void InfluenceKernel(const double *B0,
                            const double *B1,
                            const double *B2,
                            const double *B3,
                            const double *B4,
                            const double *B5,
                            const int ld3,
                            double *grid)
{  
#if defined(__SSE3__)
    // SSE3 kernel
    __m128d vB0 = _mm_loaddup_pd(B0);
    __m128d vB1 = _mm_loaddup_pd(B1);
    __m128d vB2 = _mm_loaddup_pd(B2);
    __m128d vB3 = _mm_loaddup_pd(B3);
    __m128d vB4 = _mm_loaddup_pd(B4);
    __m128d vB5 = _mm_loaddup_pd(B5);
    __m128d vg0 = _mm_load_pd(&(grid[0 * ld3]));
    __m128d vg1 = _mm_load_pd(&(grid[1 * ld3]));
    __m128d vg2 = _mm_load_pd(&(grid[2 * ld3]));          
    __m128d vv0 = _mm_mul_pd(vg0, vB0);
    vv0 = _mm_add_pd(_mm_mul_pd(vg1, vB1), vv0);          
    vv0 = _mm_add_pd(_mm_mul_pd(vg2, vB2), vv0);           
    _mm_store_pd ((double *)(&(grid[0 * ld3])), vv0);
    __m128d vv1 = _mm_mul_pd(vg0, vB1);            
    vv1 = _mm_add_pd(_mm_mul_pd(vg1, vB3), vv1);          
    vv1 = _mm_add_pd(_mm_mul_pd(vg2, vB4), vv1);  
    _mm_store_pd ((double *)(&(grid[1 * ld3])), vv1);
    __m128d vv2 = _mm_mul_pd(vg0, vB2);
    vv2 = _mm_add_pd(_mm_mul_pd(vg1, vB4), vv2);
    vv2 = _mm_add_pd(_mm_mul_pd(vg2, vB5), vv2);
    _mm_store_pd ((double *)(&(grid[2 * ld3])), vv2); 
#elif defined(__MIC__) 
    // MIC kernel
    __m512d vB0 = _mm512_mask_loadunpacklo_pd(vB0, 0xAA, B0);
    vB0 = _mm512_mask_loadunpacklo_pd(vB0, 0x55, B0);
    __m512d vB1 = _mm512_mask_loadunpacklo_pd(vB1, 0xAA, B1);
    vB1 = _mm512_mask_loadunpacklo_pd(vB1, 0x55, B1);
    __m512d vB2 = _mm512_mask_loadunpacklo_pd(vB2, 0xAA, B2);
    vB2 = _mm512_mask_loadunpacklo_pd(vB2, 0x55, B2);
    __m512d vB3 = _mm512_mask_loadunpacklo_pd(vB3, 0xAA, B3);
    vB3 = _mm512_mask_loadunpacklo_pd (vB3, 0x55, B3);
    __m512d vB4 = _mm512_mask_loadunpacklo_pd(vB4, 0xAA, B4);
    vB4 = _mm512_mask_loadunpacklo_pd(vB4, 0x55, B4);
    __m512d vB5 = _mm512_mask_loadunpacklo_pd(vB5, 0xAA, B5);
    vB5 = _mm512_mask_loadunpacklo_pd(vB5, 0x55, B5);
    __m512d vg0 = _mm512_load_pd(&(grid[0 * ld3]));
    __m512d vg1 = _mm512_load_pd(&(grid[1 * ld3]));
    __m512d vg2 = _mm512_load_pd(&(grid[2 * ld3]));            
    __m512d vv0 = _mm512_mul_pd(vg0, vB0);
    vv0 = _mm512_fmadd_pd(vg1, vB1, vv0);
    vv0 = _mm512_fmadd_pd(vg2, vB2, vv0);
    _mm512_store_pd(&(grid[0 * ld3]), vv0);
    __m512d vv1 = _mm512_mul_pd(vg0, vB1);
    vv1 = _mm512_fmadd_pd(vg1, vB3, vv1);
    vv1 = _mm512_fmadd_pd(vg2, vB4, vv1);
    _mm512_store_pd(&(grid[1 * ld3]), vv1);
    __m512d vv2 = _mm512_mul_pd(vg0, vB2);
    vv2 = _mm512_fmadd_pd(vg1, vB4, vv2);
    vv2 = _mm512_fmadd_pd(vg2, vB5, vv2);
    _mm512_store_pd(&(grid[2 * ld3]), vv2);
#else
    // scalar kernel
    double real[3];
    double imag[3];
    real[0] = grid[0 * ld3 + 0];
    imag[0] = grid[0 * ld3 + 1];
    real[1] = grid[1 * ld3 + 0];
    imag[1] = grid[1 * ld3 + 1];
    real[2] = grid[2 * ld3 + 0];
    imag[2] = grid[2 * ld3 + 1];
    grid[0 * ld3 + 0] = B0[0] * real[0] +
            B1[0] * real[1] +
            B2[0] * real[2];
    grid[0 * ld3 + 1] = B0[0] * imag[0] +
            B1[0] * imag[1] +
            B2[0] * imag[2];
    grid[1 * ld3 + 0] = B1[0] * real[0] +
            B3[0] * real[1] +
            B4[0] * real[2];
    grid[1 * ld3 + 1] = B1[0] * imag[0] +
            B3[0] * imag[1] +
            B4[0] * imag[2];
    grid[2 * ld3 + 0] = B2[0] * real[0] +
            B4[0] * real[1] +
            B5[0] * real[2];
    grid[2 * ld3 + 1] = B2[0] * imag[0] +
            B4[0] * imag[1] +
            B5[0] * imag[2];
#endif   
}


inline void InterpolateKernel(const int porder3,                                
                              const double *P,
                              const int *ind,
                              const double alpha,
                              const double *grid,
                              const int ld3,                                                            
                              const double beta,
                              double *vels)
{
    double tmp[detail::kAlignLen/sizeof(double)] __attribute__((aligned(detail::kAlignLen)));   
#if defined(__MIC__)
    __m512d vv0 = _mm512_set1_pd(0.0);
    __m512d vv1 = _mm512_set1_pd(0.0);
    __m512d vv2 = _mm512_set1_pd(0.0); 
    for (int j = 0; j < porder3; j+=detail::kSimdWidth/sizeof(double)) {
        __m512d vP = _mm512_load_pd(&(P[j]));
        __m512i vidx = _mm512_loadunpacklo_epi32(vidx, &ind[j]);          
        __m512d vg0 = _mm512_i32logather_pd(vidx, &(grid[0*ld3]), _MM_SCALE_8);
        __m512d vg1 = _mm512_i32logather_pd(vidx, &(grid[1*ld3]), _MM_SCALE_8);
        __m512d vg2 = _mm512_i32logather_pd(vidx, &(grid[2*ld3]), _MM_SCALE_8);
        vv0 = _mm512_fmadd_pd(vg0, vP, vv0);
        vv1 = _mm512_fmadd_pd(vg1, vP, vv1);
        vv2 = _mm512_fmadd_pd(vg2, vP, vv2);
    }
    tmp[0] = _mm512_reduce_add_pd(vv0);
    tmp[1] = _mm512_reduce_add_pd(vv1);
    tmp[2] = _mm512_reduce_add_pd(vv2);
    vels[0] = beta * vels[0] + alpha * tmp[0];
    vels[1] = beta * vels[1] + alpha * tmp[1];
    vels[2] = beta * vels[2] + alpha * tmp[2];
#else
    tmp[0] = 0.0;
    tmp[1] = 0.0;
    tmp[2] = 0.0;
    for (int j = 0; j < porder3; j++) {
        int idx = ind[j];
        double pvalue = P[j];
        tmp[0] += grid[0 * ld3 + idx] * pvalue;
        tmp[1] += grid[1 * ld3 + idx] * pvalue;
        tmp[2] += grid[2 * ld3 + idx] * pvalue;           
    }
    vels[0] =  beta * vels[0] + alpha * tmp[0];
    vels[1] =  beta * vels[1] + alpha * tmp[1];
    vels[2] =  beta * vels[2] + alpha * tmp[2];  
#endif
}


inline void SpreadKernel(const int porder3,
                         const double *P,
                         const int *ind,
                         const double *forces,                         
                         const int ld3,
                         double *grid)
{  
#if defined(__MIC__)
    __m512d vf0 = _mm512_set1_pd(forces[0]);
    __m512d vf1 = _mm512_set1_pd(forces[1]);
    __m512d vf2 = _mm512_set1_pd(forces[2]); 
    for (j = 0; j < porder3; j+=detail::kSimdWidth/sizeof(double)) {
        __m512d vP = _mm512_load_pd(&(P[j]));
        __m512i vidx = _mm512_loadunpacklo_epi32 (vidx, &ind[j]);        
        __m512d vg0 = _mm512_i32logather_pd(vidx, &(grid[0*ld3]), _MM_SCALE_8);
        __m512d vg1 = _mm512_i32logather_pd(vidx, &(grid[1*ld3]), _MM_SCALE_8);
        __m512d vg2 = _mm512_i32logather_pd(vidx, &(grid[2*ld3]), _MM_SCALE_8);     
        vg0 = _mm512_fmadd_pd(vP, vf0, vg0);
        vg1 = _mm512_fmadd_pd(vP, vf1, vg1);
        vg2 = _mm512_fmadd_pd(vP, vf2, vg2);
        _mm512_i32loscatter_pd(&(grid[0 * ld3]), vidx, vg0, _MM_SCALE_8);
        _mm512_i32loscatter_pd(&(grid[1 * ld3]), vidx, vg1, _MM_SCALE_8);
        _mm512_i32loscatter_pd(&(grid[2 * ld3]), vidx, vg2, _MM_SCALE_8);       
    } 
#else
    for (int j = 0; j < porder3; j++) {
        int idx = ind[j];
        double pvalue = P[j];
        grid[0 * ld3 + idx] += pvalue * forces[0];
        grid[1 * ld3 + idx] += pvalue * forces[1];
        grid[2 * ld3 + idx] += pvalue * forces[2];
    }
#endif
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


/// Creates a SPME engine
bool CreateSpmeEngine(const int npos,
                      const double box_size,
                      const double xi,
                      const int dim,
                      const int porder,
                      SpmeEngine **p_spme);

/// Destroys the SPME engine
void DestroySpmeEngine(SpmeEngine *spme);

/// Computes the reciprocal-space sum    
void ComputeSpmeRecip(const SpmeEngine *spme,
                      const int nrhs,
                      const double alpha,                              
                      const int ldin,
                      const double *vec_in,
                      const double beta,
                      const int ldout,
                      double *vec_out);

/// Updates the particles and reconstruct the P matrix   
void UpdateSpmeEngine(const double *pos, const double *rdi, SpmeEngine *spme);

} // namespace stokesdt

} // namespace detail


#endif // SPME_H_
