/**
 * @file   ewald.h
 * @brief  EwaldTable definition and operations
 */

#ifndef EWALD_H_
#define EWALD_H_


#include <cmath>


namespace stokesdt {

namespace detail {

/** 
 * @struct  EwaldTable
 * @brief   Ewald Table for efficient Ewald calculations
 */
typedef struct EwaldTable {
    /// the real-space cutoff
    double rmax;
    /// the reciprocal-space cutoff
    double kmax;    
    /// the number of real-space meshes
    int nr;    
    /// the number of reciprocal-space meshes
    int nk;
    /// the real-space <code>lx</code> table
    double *rlx_tbl  = NULL;
    /// the real-space <code>ly</code> table
    double *rly_tbl  = NULL;
    /// the real-space <code>lz</code> table
    double *rlz_tbl  = NULL;
    /// the real-space <code>ex</code> table
    double *ex_tbl   = NULL;
    /// the real-space <code>ey</code> table
    double *ey_tbl   = NULL;
    /// the real-space <code>ez</code> table
    double *ez_tbl   = NULL;
    /// the reciprocal-space <code>k</code> table
    double *k_tbl    = NULL;
    /// the reciprocal-space <code>k1</code> table
    double *k1_tbl   = NULL;
    /// the reciprocal-space <code>k2</code> table
    double *k2_tbl   = NULL;
    /// the reciprocal-space <code>k3</code> table
    double *k3_tbl   = NULL;
    /// the reciprocal-space <code>kexp</code> table
    double *kexp_tbl = NULL;
} EwaldTable;


__attribute__((vector))
inline void scalars_ewald_F(const double xi,
                            const double r,
                            const double aa,
                            const double ab,
                            double *xa,
                            double *ya)
{
    double xir = xi * r;
    double xir2 = xir * xir;
    double r2 = r * r;
    double r3 = r2 * r;
    double aa2 = aa * aa;
    double ab2 = ab * ab;
    double erfcxir = erfc(xir);
    double expxir2 = xi / sqrt(M_PI) * exp(-xir2);
    double a0 = 2.0 * (-3.0 + xir2 * 2.0) * expxir2 + erfcxir / r;
    double b0 = 2.0 * (1.0 + xir2 * -2.0) * expxir2 + erfcxir / r;
    double a2 = 4.0 / r2 * (1.0 + xir2 * (14.0 + xir2 * (-20.0 + xir2 * 4.0)))
        * expxir2 + 2.0 / r3 * erfcxir;
    double b2 = -4.0 / r2 * (3.0 + xir2 * (2.0 + xir2 * (-16.0 + xir2 * 4.0)))
        * expxir2 - 6.0 / r3 * erfcxir;
    double a10 = a0;
    double a12 = a2 * (aa2 + ab2) / 6.0;
    double a20 = b0;
    double a22 = b2 * (aa2 + ab2) / 6.0;

    // xa, ya
    *ya = 0.75 * (a10 + a12);
    *xa = 0.75 * (a20 + a22) + *ya;
}

                   
/// Constructs the mobility matrix
void EwaldKernel(const double xi, const EwaldTable *ewald_tbl,
                 const double box_size, const int npos,
                 const double *pos, const double *rdi,
                 const int ldmat, double *mat,
                 double *mat_real, double *mat_recip);

/// Computes matrix-vector product of the mobility matrix
void EwaldVectorKernel(const double xi, const EwaldTable *ewald_tbl,
                       const double box_size, const int npos,
                       const double *pos, const double *rdi,
                       const double alpha, const int ldf, const double *f,
                       const double beta, const int ldv, double *v,
                       double *v_real, double *v_recip);

/// Constructs the mobility matrix with free boundary conditions
void NonEwaldKernel(const int npos,
                    const double *pos,
                    const double *rdi,
                    const int ldm,
                    double *mat);
                    
/// Creates the Ewald tables
bool CreateEwaldTable(const double xi, const double box_size,
                      const double tol, EwaldTable **p_ewald_tbl);

/// Initializes the real-space tables
bool InitRealTable(const double xi, const double box_size,
                   const double tol, EwaldTable *tbl);

/// Initializes the reciprocal-space tables
bool InitRecipTable(const double xi, const double box_size,
                    const double tol, EwaldTable *tbl);

/// Destroys the Ewald tables
void DestroyEwaldTable(EwaldTable *ewald_tbl);


} // namespace detail

} // namespace stokesdt


#endif // EWALD_H_
