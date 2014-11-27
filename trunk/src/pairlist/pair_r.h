/**
 * @file   pair_r.h
 * @brief  PairListR class definition
 */

#ifndef VERLET_R_H_
#define VERLET_R_H_


#include <vector>
#include "pair_base.h"


namespace stokesdt{

/** 
 * @class  PairListR
 * @brief  Pair list for absolute distance
 *
 * PairListR is a tool for finding all particle pairs within
 * a given cutoff distance, which can be used to compute short-range
 * interactions between particles.
 * 
 */
class PairListR : public PairListBase {
  public:
    /**
     * @brief Class constructor
     * 
     * Constructs a new PairListR instance.
     *
     * @param[in] npos      the number of particles
     * @param[in] rdi       the array of particle radii
     * @param[in] box_size  the dimension of the simulation box
     * @param[in] cutoff    the cutoff distance
     */
    PairListR(const int npos,
                const double *rdi,
                const double box_size,
                const double cutoff);

    /// @copydoc PairListBase::~PairListR()
    virtual ~PairListR();

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
    DISALLOW_COPY_AND_ASSIGN(PairListR);

    /// @copydoc PairListBase::NumCells()
    virtual int NumCells();

    /// @copydoc PairListBase::InitSearchIdx()
    virtual void InitSearchIdx();

    /// @copydoc PairListBase::InitPairs()
    virtual int InitPairs(const double *rdi,
                          std::vector<std::vector<int> > *pairs);

    /// @copydoc PairListBase::FindPairs()
    virtual int FindPairs(const double *pos,
                          const double *rdi,
                          std::vector<std::vector<int> > *pairs);
};

} // namespace stokesdt


#endif // VERLET_LIST_H_
