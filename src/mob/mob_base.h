/**
 * @file   mob_base.h
 * @brief  the MobBase class definition
 */
 
#ifndef MOB_BASE_H_
#define MOB_BASE_H_


#include <cmath>
#include "common.h"


namespace stokesdt {

/** @class  MobBase
 *  @brief  Abstract base class for constructing mobility matrix
 */
class MobBase {
  public:
    /// Class deconstructor
    virtual ~MobBase();

    /**
     * @brief  Initializes the instance.
     * @return <code>true</code> if the instance is initialized successfully;
     *         <code>false</code> otherwise. 
     */
    virtual bool Init() = 0;

    /**
     * @brief  Updates particles.
     * 
     * Update the coordinates and radii of particles, and reconstructs the
     * mobility matrix.
     * <p>
     * This method must be called before invocating MulVector(). 
     *
     * @param[in] pos  the new coordinates of the particles
     * @param[in] rid  the radii of the particles
     */
    virtual void Update(const double *pos, const double *rdi) = 0;

    /** 
     * @brief  Multiply the mobility matrix by mutiple vectors
     *  
     * The operation is defined as
     * <pre><code>
     *     v := alpha * mob * f + beta * v
     * </code></pre> 
     *
     *  @param[in]    num_rhs  the number of the vectors to be multiplied
     *  @param[in]    alpha    the <code>alpha</code> scalar
     *  @param[in]    ldf      the leading dimension of \c f
     *  @param[in]    f        the pointer to the vectors \c f
     *  @param[in]    beta     the <code>beta</code> scalar
     *  @param[in]    ldv      the leading dimension of \c v
     *  @param[inout] v        the pointer to vectors \c v
     */
    virtual void MulVector(const int num_rhs,
                           const double alpha,
                           const int ldf,
                           const double *f,
                           const double beta,
                           const int ldv,
                           double *v) = 0;

    /// Returns the dimension of the mobility matrix
    virtual int dim() = 0;

  protected:
    /** 
     * @brief  Sole constructor. (For invocation by
     * derived class constructors, typically implicit.)
     */
    MobBase();
    
  private:
    DISALLOW_COPY_AND_ASSIGN(MobBase);
};

} // namespace stokesdt


#endif // MOB_BASE_H_
