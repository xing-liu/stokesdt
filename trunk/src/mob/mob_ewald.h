/**
 * @file   mob_ewald.h
 * @brief  the MobEwald class definition
 */

#ifndef MOB_EWALD_H_
#define MOB_EWALD_H_


#include "mob_matrix.h"
#include "ewald.h"


namespace stokesdt{

/** 
 * @class  MobEwald
 * @brief  Computes mobility matrices using Ewald summations.
 */
class MobEwald : public MobMatrix {
  public:
    /**
     * @brief  Class constructor specifying the Ewald parameter
     *
     * Constrcuts a new MobEwald instance.
     * 
     * @param[in] npos      the number of particles
     * @param[in] rdi       the radii of particles
     * @param[in] box_size  the dimension of the simulation box size
     * @param[in] tol       the requested tolerance of Ewald errors
     * @param[in] xi        the Ewald parameter
     */
    MobEwald(const int npos,
             const double *rdi,
             const double box_size,
             const double tol,
             const double xi);

    /**
     * @brief  Class constructor
     *
     * Constrcuts a new MobEwald instance.
     * 
     * @param[in] npos      the number of particles
     * @param[in] rdi       the radii of particles
     * @param[in] box_size  the dimension of the simulation box size
     * @param[in] tol       the requested tolerance of Ewald errors
     */
    MobEwald(const int npos,
             const double *rdi,
             const double box_size,
             const double tol);

    /// @copydoc MobBase::~MobBase()
    virtual ~MobEwald();

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
    DISALLOW_COPY_AND_ASSIGN(MobEwald);
    
  private:
    /// the dimension of the simulate box
    double box_size_; 
    /// the requested Ewald error
    double tol_;
    /// the Ewald parameter
    double xi_;
    /// the number of particles
    int npos_;
    /// the pointer to the dense mobility matrix
    double *mat_;
    /// the dimension of the mobility matrix
    int dim_mob_;
    /// the leading dimension of the mobility matrix
    int ldm_;
    /// the pointer to the Ewald table 
    detail::EwaldTable *ewald_tbl_;
};

} // namespace stokesdt


#endif // MOB_EWALD_H_
