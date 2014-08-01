/**
 * @file   verlet_s.h
 * @brief  VerletListS class definition
 */

#ifndef VERLET_S_H_
#define VERLET_S_H_


#include <vector>
#include "verlet_base.h"


namespace stokesdt {

/** @class  VerletListS
 *  @brief  Verlet list for normalized distance
 *
 * VerletListR is a tool for finding all particle pairs within
 * a given cutoff distance, which can be used to compute short-range
 * interactions between particles.
 * <p>
 * The normalized distance between particle
 * <code>i</code> and <code>j</code> is defined as:
 * <p>
 * <code><pre>
 *    s_ij = 2*r_ij/(a_i+a_j)
 * </pre></code>
 * where : <p>
 * <code>r_ij</code> - the distane between particle
 * <code>i</code> and <code>j</code>
 * <p>
 * <code>a_i</code>  - the radius of particle <code>i</code>
 * <p>
 * <code>a_j</code>  - the radius of particle <code>j</code>     
 */
class VerletListS : public VerletListBase {
  public:
    /**
     * @brief Class constructor
     * 
     * Constructs a new VerletListS instance.
     *
     * @param[in] npos      the number of particles
     * @param[in] rdi       the array of particle radii
     * @param[in] box_size  the dimension of the simulation box
     * @param[in] cutoff    the normalied distance cutoff
     */
    VerletListS(const int npos,
                const double *rdi,
                const double box_size,
                const double cutoff);

    /// @copydoc VerletListBase::~VerletListS()
    virtual ~VerletListS();

  private:
    /** \struct  SearchIdxS
     *  @brief   Defines the earch space for normalized distance 
     */
    struct SearchIdxS {
        /// the delta on x
        int dx;
        /// the delta on y 
        int dy;
        /// the delta on z
        int dz;
        /// the maximum radius sum of any two particles that can interact
        double max_rdi_sum;
        /// the minimun radius sum of any two particles that can interact
        double min_rdi_sum;
    };

    /// the search indices
    std::vector<SearchIdxS> search_idx_;
      
  private:
    DISALLOW_COPY_AND_ASSIGN(VerletListS);

    /// @copydoc VerletListBase::NumCells()
    virtual int NumCells();

    /// @copydoc VerletListBase::InitSearchIdx()
    virtual void InitSearchIdx();

    /// @copydoc VerletListBase::InitPairs()
    virtual int InitPairs(const double *rdi,
                          std::vector<std::vector<int> > *pairs);

    /// @copydoc VerletListBase::FindPairs()
    virtual int FindPairs(const double *pos,
                          const double *rdi,
                          std::vector<std::vector<int> > *pairs);
};

} // namespace stokesdt


#endif // VERLET_S_H_