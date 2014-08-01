/**
 * @file   brwn_base.h
 * @brief  BrwnBase class definition
 */

#ifndef BRWN_BASE_H_
#define BRWN_BASE_H_


#include "common.h"
#include "mob_base.h"

 
namespace stokesdt {

/**
 * @class  BrwnBase
 * @brief  Abstract base class for computing brownian displacement vectors
 */
class BrwnBase {
  public:
    /// @brief  Class deconstructor
    virtual ~BrwnBase();

    /**
     * @brief  Initializes the instance.
     * @return <code>true</code> if the instance is initialized successfully;
     *         <code>false</code> otherwise. 
     */
    virtual bool Init() = 0;

    /** 
     * @brief  Computes the Brownian displacement vectors
     *         for a given MobBase matrix and mutiple random vectors.
     *
     * @param[in]  mob      the pointer to the MobBase instance
     * @param[in]  num_rhs  the number of random vetors to be computed
     * @param[in]  ldz      the leading dimension of \c z
     * @param[in]  z        the pointer to the multiple random vectors
     * @param[in]  ldy      the leading dimension of of \c y
     * @param[out] y        the pointer to the Brownian displacements to
     *                      be computed
     */
    virtual void Compute(MobBase *mob, const int num_rhs,
                         const int ldz, const double *z,
                         const int ldy, double *y) = 0;
    
  protected:
    /** 
     * @brief  Sole constructor. (For invocation by
     * derived class constructors, typically implicit.)
     */
    BrwnBase();
    
  private:
    DISALLOW_COPY_AND_ASSIGN(BrwnBase);
};

} // namespace stokesdt


#endif // BRWN_BASE_H_
