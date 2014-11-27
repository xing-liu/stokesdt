/**
 * @file   force_steric.h
 * @brief  StericForce class definition
 */
 
#ifndef FORCE_STERIC_H_
#define FORCE_STERIC_H_


#include <vector>
#include "pair_s.h"
#include "force_base.h"


namespace stokesdt {

/** @class  StericForce
 *  @brief  Computes steric forces using PairListS
 */
class StericForce : public ForceBase {
  public:
    /**
     * @brief Class constructor
     *
     * Constructs a new StericForce instance.
     * <p>
     * PairListS is used to compute the short-range steric forces. 
     * <p>
     * The periodic bondary condition is assumed.
     * The following harmonic function is
     * used to model the steric force between particis
     * <code>i</code> and <code>j</code>:
     * <p>
     * <code><pre>
     *    f_ij = k0*(2*r_ij/(a_i+a_j)-r0), if 2*r_ij/(a_i+a_j)-r0 <= 0;
     *         = 0,                            otherwise.
     * </pre></code>
     * where : <p>
     * <code>k0</code>   - the steric force constant
     * <p>
     * <code>r_ij</code> - the distane between particle
     * <code>i</code> and <code>j</code>
     * <p>
     * <code>a_i</code>  - the radius of particle <code>i</code>
     * <p>
     * <code>a_j</code>  - the radius of particle <code>j</code>
     * <p>
     * <code>r0</code>   - the steric cutoff
     *
     * @param[in] npos       the number of particles
     * @param[in] rdi        th array of particle radii
     * @param[in] box_size   the dimension of the simulation box
     * @param[in] steric_r0  the cutoff of the steric force
     * @param[in] steric_k0  the steric force constant
     */
    StericForce(const int npos,
                const double *rdi,
                const double box_size, 
                const double steric_r0,
                const double steric_k0);

    /// @copydoc ForceBase::~ForceBase()
    virtual ~StericForce();

    /// @copydoc ForceBase::Init()
    virtual bool Init();

    /// @copydoc ForceBase::Accumulate()
    virtual void Accumulate(const double *pos, const double *rdi,
                            double *f);
    
  private:
    DISALLOW_COPY_AND_ASSIGN(StericForce);

    /// the kernel for computing steric forces
    void ForceKernel(const double *pos,
                     const double *rdi,
                     double *f);

  private:
    /// the steric cutoff
    double steric_r0_;
    /// the steric coefficient
    double steric_k0_;
    /// the dimension of the simulation box
    double box_size_;
    /// the PairListS instance for computing steric forces
    PairListS pair_list_;
    /// the array of row pointers for storing the interactions
    std::vector<int> rowptr_;
    /// the array of column indices for storing the interations
    std::vector<int> colidx_;
};

} // namespace stokesdt


#endif // FORCE_STERIC_H_