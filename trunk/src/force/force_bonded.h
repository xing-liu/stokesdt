/**
 * @file   force_bonded.h
 * @brief  BondedForce class definition
 */
 
#ifndef FORCE_BOND_H_
#define FORCE_BOND_H_


#include <vector>
#include <tuple>
#include "force_base.h"


namespace stokesdt {

/** @class  BondedForce
 *  @brief  Computes bond forces
 */
class BondedForce : public ForceBase {
  public:
    /**
     * @brief Class constructor
     *
     * Constructs a new BondedForce instance.
     * <p>
     * Each bond can be represented by a 4-tuple:
     * <p>
     * <code><pre>
     *     <id_a, id_b, r0, k0>
     * </pre></code>
     * where:
     * <p>
     * <code>id_a</code> and <code>id_b</code> are the particle indices.
     * <p>
     * <code>r0</code> is the equilibrium bond length.
     * <p>
     * <code>k0</code> is the bond force constant.
     * <p>
     * The following harmonic function is
     * used to he bond force between particis <code>i</code> and <code>j</code>:
     * <p>
     * <code><pre>
     *    f_ij = k0_ij*(r_ij-r0_ij)
     * </pre></code>
     * where:
     * <p>
     * <code>r_ij</code> is the distane between
     * particle <code>i</code> and <code>j</code>
     * <p>
     * <code>k0_ij</code> is the bond force constant between
     * particle <code>i</code> and <code>j</code>
     * <p>
     * <code>r0_ij</code> is the equilibrium bond length between
     * particle <code>i</code> and <code>j</code>
     *
     * @param[in] npos       the number of particles
     * @param[in] num_bonds  the number of bonds to be computed
     * @param[in] bond_id_a  the array of \c id_a
     * @param[in] bond_id_b  the array of \c id_b
     * @param[in] bond_r0    the array of \c r0
     * @param[in] bond_k0    the array of \c k0
     */
    BondedForce(const int npos,
                const int num_bonds,
                const int *bond_id_a,
                const int *bond_id_b,
                const double *bond_r0,
                const double *bond_k0);

    /// @copydoc ForceBase::~ForceBase()
    virtual ~BondedForce();

    /// @copydoc ForceBase::Init()
    virtual bool Init();

    /// @copydoc ForceBase::Accumulate()
    virtual void Accumulate(const double *pos, const double *rdi,
                            double *f);
    
  private:
    DISALLOW_COPY_AND_ASSIGN(BondedForce);

  private:
    /// the array of the bond parameters
    std::vector<std::tuple<int,int,double,double> > bonds_;
    /// the array of the bond row pointers
    std::vector<int> rowptr_;
    /// the array of the bond column indices 
    std::vector<int> colidx_;
    /// the array of the equilibrium bond lengths
    std::vector<double> bond_r0_;
    /// the array of the bond force constants
    std::vector<double> bond_k0_;
};


namespace detail {
    
inline bool BondCompare1(const std::tuple<int,int,double,double> &lhs,
                         const std::tuple<int,int,double,double> &rhs)
{
  return (std::get<1>(lhs) < std::get<1>(rhs));
}

inline bool BondCompare2(const std::tuple<int,int,double,double> &lhs,
                         const std::tuple<int,int,double,double> &rhs)
{
    return (std::get<2>(lhs) < std::get<2>(rhs));
}

} // namespace detail

} // namespace stokesdt

#endif // FORCE_BOND_H_