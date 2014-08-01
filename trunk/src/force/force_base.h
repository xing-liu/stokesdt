/**
 * @file   force_base.h
 * @brief  ForceBase class definition
 */

#ifndef FORCE_BASE_H_
#define FORCE_BASE_H_


#include "common.h"


namespace stokesdt {

/** 
 * @class  ForceBase
 * @brief  Abstract base class for computing forces
 */
class ForceBase {
  public:
    /** 
     * @brief  Class constructor specifying the number of
     *         particles to be computed
     *
     * @param[in] npos  the number of particles  
     */
    ForceBase(int npos);
        
    /// Class deconstructor
    virtual ~ForceBase();

    /**
     * @brief  Initializes the instance.
     * @return <code>true</code> if the instance is initialized successfully;
     *         <code>false</code> otherwise. 
     */
    virtual bool Init() = 0;

    /**
     * @brief  Computes the force vector for a given position vector
     *         and a radius vector.
     * 
     * @param[in]  pos  the position vector
     * @param[in]  rdi  the radius vector
     * @param[out] f    the pointer to the force vector
     */
    void Compute(const double *pos, const double *rdi, double *f);

    /**
     * @brief  Computes forces and accumulate them into a given force vector.
     * 
     * This function is equalavent to the following. <p>
     * <code><pre>
     *     Compute(pos, rdi, f0); <p>
     *     f += f0;
     * </pre></code>
     * 
     * @see Compute()
     *
     * @param[in]    pos  the position vector
     * @param[in]    rdi  the radius vector
     * @param[inout] f    the pointer to the force vector
     */
    virtual void Accumulate(const double *pos, const double *rdi,
                            double *f) = 0;

  protected:
    /**
     * @brief  the access function returning the value of <code>npos_</code>
     *         (For invocation by derived classes only.)
     */
    int get_npos() {return npos_;};
            
  private:
    DISALLOW_COPY_AND_ASSIGN(ForceBase);

  private:
    /// the number of particles
    int npos_;
};

} // namespace stokesdt


#endif // FORCE_BASE_H_
