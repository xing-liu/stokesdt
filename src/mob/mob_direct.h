/**
 * @file   mob_direct.h
 * @brief  the MobDirect class definition
 */

#ifndef MOB_DIRECT_H_
#define MOB_DIRECT_H_


#include "mob_matrix.h"
#include "ewald.h"


namespace stokesdt{

/** 
 * @class  MobDirect
 * @brief  Computes mobility matrices with free boundary conditions.
 */
class MobDirect : public MobMatrix {
  public:
    /**
     * @brief  Class constructor specifying the Ewald parameter
     *
     * Constrcuts a new MobDirect instance.
     * 
     * @param[in] npos      the number of particles
     * @param[in] rdi       the radii of particles
     */
    MobDirect(const int npos,
              const double *rdi);

    /// @copydoc MobBase::~MobBase()
    virtual ~MobDirect();

    /// @copydoc MobBase::Init()
    virtual bool Init();

    /// @copydoc MobBase::Update()
    virtual void Update(const double *pos, const double *rdi);

    /// @copydoc MobBase::MulVector()
    virtual void MulVector(const int num_rhs,
                           const double alpha,
                           const int ldvec_in,
                           const double *vec_in,
                           const double beta,
                           const int ldvec_out,
                           double * vec_out);

    /// @copydoc MobMatrix::GetMatrix()
    virtual void GetMatrix(const int ldm, double *mat);

    /// @copydoc MobBase::dim()
    virtual int dim() {return dim_mob_;};
    
  private:
    DISALLOW_COPY_AND_ASSIGN(MobDirect);
    
  private:
    /// the number of particles
    int npos_;
    /// the pointer to the dense mobility matrix
    double *mat_;
    /// the dimension of the mobility matrix
    int dim_mob_;
    /// the leading dimension of the mobility matrix
    int ldm_;
};

} // namespace stokesdt


#endif // MOB_DIRECT_H_
