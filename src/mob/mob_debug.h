/**
 * @file   mob_debug.h
 * @brief  the MobDebug class definition
 */

#ifndef MOB_DEBUG_H_
#define MOB_DEBUG_H_


#include "ewald.h"
#include "common.h"


namespace stokesdt{

/** @class  MobDebug
 *  @brief  Computes mobility matrix using Ewald summations (for debug).
 */
class MobDebug {
  public:
    typedef enum MobDebugType_{
        EWALD       = 0,
        EWALD_REAL  = 1,
        EWALD_RECIP = 2
    } MobDebugType;
    
    /// Class constructor specifying the Ewald parameter
    MobDebug(const int npos,
             const double *rdi,
             const double box_size,
             const double tol,
             const double xi,
             MobDebugType mode = EWALD);

    /// Class Constructor
    MobDebug(const int npos,
             const double *rdi,
             const double box_size,
             const double tol,
             MobDebugType mode = EWALD);

    /// Class deconstructor
    virtual ~MobDebug();

    /// Initializes the instance.
    bool Init();

    /// Returns the dimension of the mobility matrix.
    int dim() {return dim_mob_;};

    void MulVector(const double *pos,
                   const double *rdi,
                   const int num_rhs,
                   const double alpha,
                   const int ldf,
                   const double *f,
                   const double beta,
                   const int ldv,
                   double *v);

    /// Multiplies the real-space part of the mobility matrix.
    void RealMulVector(const double *pos,
                       const double *rdi,
                       const int num_rhs,
                       const double alpha,
                       const int ldf,
                       const double *f,
                       const double beta,
                       const int ldv,
                       double *v);

    /// Multiplies the reciprocal-space part of the mobility matrix.
    void RecipMulVector(const double *pos,
                        const double *rdi,
                        const int num_rhs,
                        const double alpha,
                        const int ldf,
                        const double *f,
                        const double beta,
                        const int ldv,
                        double *v);

    /// Returns the mobility matrix into a local buffer
    void GetMatrix(const double *pos, const double *rdi,
                   const int ldm, double *mat);
    
  private:
    DISALLOW_COPY_AND_ASSIGN(MobDebug);
    
  private:
    /// the Ewald mode, REAL
    MobDebugType mode_;
    /// the dimension of the simulate box
    double box_size_; 
    /// the requested tolerance for the Ewald errors
    double tol_;
    /// the Ewald parameter
    double xi_;
    /// the number of particles
    int npos_;
    /// the dimension of the mobility matrix
    int dim_mob_ = 0;
    /// the Ewald tables
    detail::EwaldTable *ewald_tbl_ = NULL;
};

} // namespace stokesdt


#endif // MOB_DEBUG_H_
