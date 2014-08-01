/**
 * @file   mob_matrix.h
 * @brief  the MobMatrix class definition
 */

#ifndef MOB_MATRIX_H_
#define MOB_MATRIX_H_


#include "mob_base.h"


namespace stokesdt{

/** @class  MobMatrix
 *  @brief  Abstract base class for explicitly constructing
 *          and storing mobilit matrices
 */
class MobMatrix : public MobBase {
  public:
    /// @copydoc MobBase::~MobBase()
    virtual ~MobMatrix();

    /// @copydoc MobBase::Init()
    virtual bool Init() = 0;

    /// @copydoc MobBase::Update()
    virtual void Update(const double *pos, const double *rdi) = 0;

    /// @copydoc MobBase::MulVector()
    virtual void MulVector(const int num_rhs,
                           const double alpha,
                           const int ldvec_in,
                           const double * vec_in,
                           const double beta,
                           const int ldvec_out,
                           double * vec_out) = 0;

    /** @brief  Returns the dense mobility matrix into a local buffer
     *
     *  @param[in]  ldm  the leading dimension of the local buffer
     *  @param[out] mat  the pointer to the local buffer where the data goes
     */
    virtual void GetMatrix(const int ldm, double *mat) = 0;

    /// @copydoc MobBase::dim()
    virtual int dim() = 0;

  protected:
    /// @copydoc MobBase::MobBase()
    MobMatrix();
        
  private:
    DISALLOW_COPY_AND_ASSIGN(MobMatrix);
};

} // namespace stokesdt


#endif // MOB_MATRIX_H_
