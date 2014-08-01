/**
 * @file   brwn_chol.h
 * @brief  BrwnChol class definition
 */
 
#ifndef BRWN_CHOL_H_
#define BRWN_CHOL_H_


#include "brwn_base.h"


namespace stokesdt {

/** 
 *  @class  BrwnChol
 *  @brief  Computes Brownian displacement vectors using Cholesky factorization.
 */
class BrwnChol : public BrwnBase {
  public:
    /** 
     * @brief  Class constructor
     *
     * Constructs a new BrwnChol instance that uses Cholesky factorization 
     * to compute Brownian displacement vectors for a given MobBase instance.
     * <p>
     * The <code>dim</code> argument specifies the dimension of
     * the mobility matrix represented by the MobBase instance.
     *
     * @param[in] dim  the dimesnion of the MobBase instance to be computed
     */
    BrwnChol(const int dim);

    /// \copydoc  BrwnBase::~BrwnBase()
    virtual ~BrwnChol();

    /// \copydoc  BrwnBase::Init()
    virtual bool Init();

    /// \copydoc  BrwnBase::Compute()
    virtual void Compute(MobBase *mob, const int num_rhs,
                         const int ldz, const double *z,
                         const int ldy, double *y);

  private:
    DISALLOW_COPY_AND_ASSIGN(BrwnChol);
    
  private:
    /// the dimension of the mobility matrix
    int dim_;
    /// the leading dimension of the mobility matrix
    int ldm_;
    /// the Cholesky factor <code>L</code> of the mobility matrix
    double *cholmat_;
};

} // namespace stokesdt


#endif // BRWN_CHOL_H_
