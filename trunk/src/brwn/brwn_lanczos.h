/**
 * @file   brwn_lanczos.h
 * @brief  BrwnLanczos class definition
 */

#ifndef BRWN_LANCZOS_H_
#define BRWN_LANCZOS_H_


#include "brwn_base.h"


namespace stokesdt {

/** 
 * @class BrwnLanczos
 * @brief Compute Brownian displacement vectors using
 *        the Krylov subspace method.
 */
class BrwnLanczos : public BrwnBase {
  public:
    /** 
     * @brief  Class constructor
     *
     * Constructs a new BrwnLanczos instance that uses the Krylov subspace 
     * method (the Lanczos method) to compute Brownian displacement vectors
     * for a given MobBase instance.
     * <p>
     * The <code>dim</code> argument specifies the dimension of
     * the mobility matrix represented by the MobBase instance.
     * <p>
     * The Brownian displacement vectors are divided into groups, and
     * the block Lanczos method is used to compute all vectors of a group
     * simultaneously. The argument <code>max_nrhs</code>
     * specifies the maximum size of a group, that is the maximum number
     * of vectors that can be computed together.
     *
     * @param[in] dim        the dimesnion of the MobBase instance
     *                       to be computed
     * @param[in] max_iters  the maximum number of iterations.
     * @param[in] max_nrhs   the maximum number of Brownian displacement
     *                       vectors that can be computed simultaneously
     * @param[in] tol        the requested tolerance
     */
    BrwnLanczos(const int dim, const int max_iters,
                const int max_nrhs, const double tol);

    /// @copydoc BrwnBase::~BrwnBase()
    virtual ~BrwnLanczos();

    /// @copydoc  BrwnBase::Init()
    virtual bool Init();

    /// @copydoc  BrwnBase::Compute()
    virtual void Compute(MobBase *mob, const int nrhs,
                         const int ldy, const double *y,
                         const int ldz, double *z);

  private:
    DISALLOW_COPY_AND_ASSIGN(BrwnLanczos);

    /**
     * /brief  Computes Lanczos for a single vector
     * @param[in]  mob  the pointer to the MobBase instance
     * @param[in]  z    the pointer to the random vector
     * @param[out] y    the pointer to the Brownian displacement vector to
     *                  be computed
     * @return          the iteration number at which the Lanczos iteration
     *                  terminates
     */
    int Lanczos(MobBase *mob, const double *z, double *y);

    /**
     * /brief  Computes block Lanczos for multiple vector
     * @param[in]  mob      the pointer to the MobBase instance
     * @param[in]  num_rhs  the number of random vetors to be computed
     * @param[in]  ldz      the leading dimension of \c z
     * @param[in]  z        the pointer to the multiple random vectors
     * @param[in]  ldy      the leading dimension of of \c y
     * @param[out] y        the pointer to the Brownian displacements to
     *                      be computed
     * @return              the iteration number at which the Lanczos iteration
     *                      terminates
     */
    int BlockLanczos(MobBase *mob,
                     const int num_rhs,
                     const int ldz,
                     const double *z,
                     const int ldy,
                     double *y);
                                                       
  private:
    /// the maximum number of iterations
    int max_iters_;
    /// the maximum number of right hand sides
    int max_nrhs_;
    /// the request tolerance
    double tol_;
    /// the dimension of the mobility matrix
    int dim_;
    /// the buffer storing the Krolev sequence
    double *v_;
    /// the leading dimension of \c v_
    int ldv_;
    /// the buffer for old y
    double *y_old_;
};

} // namespace stokesdt


#endif // BRWN_LANCZOS_H_
