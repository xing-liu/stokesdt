/**
 * @file   verlet_r.h
 * @brief  VerletListR class definition
 */

#ifndef VERLET_R_H_
#define VERLET_R_H_


#include <vector>
#include "verlet_base.h"


namespace stokesdt{

/** 
 * @class  VerletListR
 * @brief  Verlet list for absolute distance
 *
 * VerletListR is a tool for finding all particle pairs within
 * a given cutoff distance, which can be used to compute short-range
 * interactions between particles.
 * 
 */
class VerletListR : public VerletListBase {
  public:
    /**
     * @brief Class constructor
     * 
     * Constructs a new VerletListR instance.
     *
     * @param[in] npos      the number of particles
     * @param[in] rdi       the array of particle radii
     * @param[in] box_size  the dimension of the simulation box
     * @param[in] cutoff    the cutoff distance
     */
    VerletListR(const int npos,
                const double *rdi,
                const double box_size,
                const double cutoff);

    /// @copydoc VerletListBase::~VerletListR()
    virtual ~VerletListR();

  private:
    /** \struct  SearchIdxR
     *  @brief   Defines the search space for absolute distance 
     */
    struct SearchIdxR {
        /// the delta on x
        int dx;
        /// the delta on y 
        int dy;
        /// the delta on z
        int dz;
        /// include all?
        bool all;
    };

    /// the search indices
    std::vector<SearchIdxR> search_idx_;
        
  private:
    DISALLOW_COPY_AND_ASSIGN(VerletListR);

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


#endif // VERLET_LIST_H_
