/**
 * @file   mob_spme.h
 * @brief  MobSpme class definition
 */

#ifndef MOB_SPME_H_
#define MOB_SPME_H_


#include "mob_base.h"
#include "spme.h"
#include "sparse.h"
#include "pair_r.h"


namespace stokesdt {

/** @class  MobSpme
 *  @brief  Computes mobility matrix using SPME.
 */
class MobSpme : public MobBase {
  public:
    /** 
     * @brief Class constructor
     * 
     */
    MobSpme(const int npos,
            const double *rdi,
            const double box_size,
            const double xi,
            const double rmax,
            const int dim,
            const int porder);

    /// @copydoc MobBase::~MobBase()
    virtual ~MobSpme();

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

    /** @brief  Multiply the real-space part of the mobility matrix
 *          by a block of vectors. The operation is defined as
 *            vec_out := alpha*mob*vec_in + beta*vec_out,   
 *
 *  @param[in]  num_rhs    Number of the vectors
 *  @param[in]  alpha      Constant alpha
 *  @param[in]  ldvec_in   Leading dimension of vec_in
 *  @param[in]  vec_in     Vector input
 *  @param[in]  beta       Constant beta
 *  @param[in]  ldvec_out  Leading dimension of vec_out
 *  @param[out] vec_out    Vector output
 */
    void RealMulVector(const int num_rhs,
                       const double alpha,
                       const int ldvec_in,
                       const double *vec_in,
                       const double beta,
                       const int ldvec_out,
                       double * vec_out);

    /** @brief  Multiply the real-space part of the mobility matrix
 *          by a block of vectors. The operation is defined as
 *            vec_out := alpha*mob*vec_in + beta*vec_out,   
 *
 *  @param[in]  num_rhs    Number of the vectors
 *  @param[in]  alpha      Constant alpha
 *  @param[in]  ldvec_in   Leading dimension of vec_in
 *  @param[in]  vec_in     Vector input
 *  @param[in]  beta       Constant beta
 *  @param[in]  ldvec_out  Leading dimension of vec_out
 *  @param[out] vec_out    Vector output
 */
    void RecipMulVector(const int num_rhs,
                        const double alpha,
                        const int ldvec_in,
                        const double *vec_in,
                        const double beta,
                        const int ldvec_out,
                        double * vec_out);

    /// @copydoc MobBase::dim()
    virtual int dim() {return dim_mob_;};
                       
  private:
    DISALLOW_COPY_AND_ASSIGN(MobSpme);

    /// Computes the real-space sums
    void BuildSparseReal(const double *pos, const double *rdi);

  private:
    /// the number of particles
    int npos_;
    /// the dimension of the mobility matrix 
    int dim_mob_;
    /// the array of particle radii
    std::vector<double> rdi_;
    /// the dimension of the simulation box
    double box_size_;
    /// the dimension of the FFT mesh
    int dim_;
    /// the interpolation order
    int porder_;
    /// the Ewald parameter
    double xi_;
    /// the real-space cutoff
    double rmax_;
    /// the pointer to the real-space mobility matrix
    detail::SparseMatrix *real_mat_;
    /// the pointer to the SPME engine
    detail::SpmeEngine *spme_;
    /// the PairListR instance for computing real-space sums
    PairListR pair_list_;
};

} // namespace stokesdt


#endif // MOB_SPME_H_
