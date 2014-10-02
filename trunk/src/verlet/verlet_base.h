/**
 * @file   verlet_base.h
 * @brief  VerletListBase class definition
 */

#ifndef VERLET_BASE_H_
#define VERLET_BASE_H_


#include <vector>
#include "common.h"


namespace stokesdt {

/**
 * @class  VerletListBase
 * @brief  Abstract base class for Verlet list
 */
class VerletListBase {
  public:
    /// Class deconstructor
    virtual ~VerletListBase();

    /**
     * @brief  Initializes the instance.
     * @return the estimated total number of particle pairs
     *         with the given cutoff distance; 0: failed to initialize
     */
    int Init();
    
    /**
     * @brief  Constructs the cell list.
     *
     * Update the particle positions and re-builds the cell list.
     * <p>
     * This method must be called before invocating GetInteractions().
     *
     * @param[in] pos  the array of particle coordinates
     * @return         the total number of particle pairs
     *                 with the given cutoff distance 
     */
    int Build(const double *pos);

    /**
     * @brief  Returns the particle pairs ith the given cutoff distance. 
     *
     * The particle pairs are stored in the CSR format,
     * which can be represented by two arrays:
     * <code>rowptr</code> and <code>colidx</code>. 
     *
     * @param[out] rowptr  the array of row pointers
     * @param[out] colidx  the array of column indices
     */
    void GetPairs(int *rowptr, int *colidx);

  protected:
    /**
     * @brief  Class constructor (For invocation by derived classes only.)
     * 
     * @param[in] npos      the number of particles
     * @param[in] rdi       the array of particle radii
     * @param[in] box_size  the dimension of the simulation box
     * @param[in] cutoff    the cutoff distance
     */
    VerletListBase(const int npos,
                   const double *rdi,
                   const double box_size,
                   const double cutoff);

    /**
     * @brief  the access function returning the value of <code>head_</code>
     *         (For invocation by derived classes only.)
     */
    const std::vector<int> &get_head() {return head_;};

    /**
     * @brief  the access function returning the value of <code>next_</code>
     *         (For invocation by derived classes only.)
     */
    const std::vector<int> &get_next() {return next_;};

    /**
     * @brief  the access function returning the value of <code>cidx_</code>
     *         (For invocation by derived classes only.)
     */
    const std::vector<size_t> &get_cidx() {return cidx_;};

    /**
     * @brief  the access function returning the value of <code>pos0_</code>
     *         (For invocation by derived classes only.)
     */
    const std::vector<double> &get_pos0() {return pos0_;};

    /**
     * @brief  the access function returning the value of <code>rdi_</code>
     *         (For invocation by derived classes only.)
     */
    const std::vector<double> &get_rdi() {return rdi_;};

    /**
     * @brief  the access function returning the value of <code>nc1_</code>
     *         (For invocation by derived classes only.)
     */
    int get_nc1() {return nc1_;};

    /**
     * @brief  the access function returning the value of <code>box_size_</code>
     *         (For invocation by derived classes only.)
     */
    double get_box_size() {return box_size_;};

    /**
     * @brief  the access function returning the value
     *         of <code>cell_size_</code>
     *         (For invocation by derived classes only.)
     */
    double get_cell_size() {return cell_size_;};

    /**
     * @brief  the access function returning the value of <code>cutoff_</code>
     *         (For invocation by derived classes only.)
     */
    double get_cutoff() {return cutoff_;};

    /**
     * @brief  the access function returning the value of <code>npos_</code>
     *         (For invocation by derived classes only.)
     */
    int get_npos() {return npos_;};

    /**
     * @brief  the access function returning the value
     *         of <code>max_radius_</code>
     *         (For invocation by derived classes only.)
     */
    double get_max_radius() {return max_radius_;};

    /**
     * @brief  the access function returning the value
     *         of <code>min_radius_</code>
     *         (For invocation by derived classes only.)
     */
    double get_min_radius() {return min_radius_;};
        
  private:
    DISALLOW_COPY_AND_ASSIGN(VerletListBase);

    /**
     * @brief  Constructs the cell list
     *
     * @param[in] pos  the array of particle coordiantes
     * @param[in] rdi  the array of particle radii
     */
    void ConstructCellList(const double *pos, const double *rdi);

    /**
     * @brief  Finds all particle pairs within the given cut-off distance
     *
     * @param[in]  pos    the array of particle coordiantes
     * @param[in]  rdi    the array of particle radii
     * @param[out] pairs  the array of particle pairs
     */
    virtual int FindPairs(const double *pos,
                          const double *rdi,
                          std::vector<std::vector<int> > *pairs) = 0;

    /// the access function returning the value of <code>nc1_</code>
    virtual int NumCells() = 0;

    /// Initializes the search space
    virtual void InitSearchIdx() = 0;

    /**
     * @brief  Initializes the array of particle pairs
     *
     * @param[in]  rdi       the array of particle radii
     * @param[out] interact  the array of interactions
     *
     * @return               the estimated total number of particle pairs
     *                       within the given cutoff
     */
    virtual int InitPairs(const double *rdi,
                          std::vector<std::vector<int> > *pairs) = 0;
    
  private:
    /// the number of particles
    int npos_;
    /// the array of particle radii
    std::vector<double> rdi_;
    /// the dimension of the simulation box size
    double box_size_;
    /// the Verletlist cutoff
    double cutoff_;
    /// the number of cells in 1D
    int nc1_;
    /// the dimension of cells
    double cell_size_;
    /// the cell list, head
    std::vector<int> head_;
    /// the cell list, next
    std::vector<int> next_;
    /// the cell list, cidx
    std::vector<size_t> cidx_;
    /// the array of normalized positions
    std::vector<double> pos0_;
    /// the maximum radius of particles
    double max_radius_;
    /// the minimum radius of particles
    double min_radius_;
    /// the array of particle pairs
    std::vector<std::vector<int> > pairs_;
};

} // namespace stokesdt


#endif // VERLET_BASE_H_