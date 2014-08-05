/**
 * @file   force_bonded.cc
 * @brief  BondedForce class implementation
 */

#include <cmath>
#include <algorithm>
#include "force_bonded.h"
#include "log.h"
#include "profile.h"


namespace stokesdt {

BondedForce::BondedForce(const int npos,
                         const int num_bonds,
                         const int *bond_id_a,
                         const int *bond_id_b,
                         const double *bond_r0,
                         const double *bond_k0)
    : ForceBase(npos)                    
{
    std::tuple<int,int,double,double> bond_info;
    for (int i = 0; i < num_bonds; i++) {
        bond_info = std::make_tuple(bond_id_a[i], bond_id_b[i],
                                    bond_r0[i], bond_k0[i]);
        bonds_.push_back(bond_info);
    }
}


BondedForce::~BondedForce()
{

}


bool BondedForce::Init()
{
    // check inputs and assign rowptr_
    int npos = get_npos();
    rowptr_.assign (npos + 1, 0);
    int num_bonds = bonds_.size();
    LOG(3, "\n        Initializes BondedForce\n");
    LOG(3, "        -----------------------\n");    
    LOG(3, "number-bonds = %d\n", num_bonds);
    for (int i = 0; i < num_bonds; i++) {
        int id_a = std::get<0>(bonds_[i]);
        int id_b = std::get<1>(bonds_[i]);
        double r0 = std::get<2>(bonds_[i]);
        double k0 = std::get<3>(bonds_[i]);
        if (id_a < 0 || id_a >= npos ||
            id_b < 0 || id_b >= npos ||
            r0 <= 0.0 || k0 <= 0.0) {
            LOG_ERROR("The bond is invalid: <%d, %d, %g, %g>\n",
                id_a, id_b, r0, k0);
            return false;
        } else {
            rowptr_[id_a + 1]++;
            rowptr_[id_b + 1]++;
        }
    }    

    // reduce rowptr_
    for (int i = 0; i < npos; i++) {
        rowptr_[i + 1] += rowptr_[i];
    }
    
    // sort bonds
    std::sort(bonds_.begin(), bonds_.end(), detail::BondCompare1);
    for (int i = 0; i < npos; i++) {
        std::sort(bonds_.begin() + rowptr_[i],
                  bonds_.begin() + rowptr_[i + 1], detail::BondCompare2);
    }
       
    // assign arrays
    colidx_.reserve(num_bonds);
    bond_r0_.reserve(num_bonds);
    bond_k0_.reserve(num_bonds);
    for (int i = 0; i < num_bonds; i++) {
        int id_b = std::get<1>(bonds_[i]);
        double r0 = std::get<2>(bonds_[i]);
        double k0 = std::get<3>(bonds_[i]);
        colidx_.push_back(id_b);
        bond_r0_.push_back(r0);
        bond_k0_.push_back(k0);        
    }

    bonds_.clear();

    return true;
}


void BondedForce::Accumulate(const double *pos, const double *rdi, double *f)
{
    START_TIMER(detail::BONDED_TICKS);
    
    int npos = get_npos();
    #pragma omp parallel for
    for (int i = 0; i < npos; i++) {
        const double x1 = pos[3 * i + 0];
        const double y1 = pos[3 * i + 1];
        const double z1 = pos[3 * i + 2];
        const int start = rowptr_[i];
        const int end = rowptr_[i + 1];
        for (int k = start; k < end; k++) {
            const int j = colidx_[k];
            if (i == j)
                continue;
            const double r0 = bond_r0_[k];
            const double kbond = bond_k0_[k];
            // compute force
            const double x2 = pos[3 * j + 0];
            const double y2 = pos[3 * j + 1];
            const double z2 = pos[3 * j + 2];
            // compute distance
            const double rx = x1 - x2;
            const double ry = y1 - y2;
            const double rz = z1 - z2;            
            const double rr = rx * rx + ry * ry + rz * rz;
            const double r = sqrt(rr);
            const double ex = rx / r;
            const double ey = ry / r;
            const double ez = rz / r;
            const double force1 = kbond * (r0 - r);
            f[3 * i + 0] +=  force1 * ex;
            f[3 * i + 1] +=  force1 * ey;
            f[3 * i + 2] +=  force1 * ez;
        }
    }

    STOP_TIMER(detail::BONDED_TICKS);
}

} //namespace stokesdt
