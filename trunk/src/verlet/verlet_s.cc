/**
 * @file   verlet_s.cc
 * @brief  VerletListS class implementation
 */

#include <algorithm>
#include <cmath>
#include <omp.h>
#include <float.h>
#include "verlet_s.h"
#include "common.h"


namespace stokesdt {

int VerletListS::NumCells()
{
    double box_size = get_box_size();
    double cutoff = get_cutoff();
    double max_radius = get_max_radius();
    int nc1 = (int)(box_size/(cutoff*max_radius/2.0));
    nc1 = (nc1 == 0 ? 1 : nc1);
    return nc1;
}
    

void VerletListS::InitSearchIdx()
{
    int nc1 = get_nc1();
    double cell_size = get_cell_size();
    double cutoff = get_cutoff();
    double max_radius_ = get_max_radius();
    size_t nc2 = size_t(nc1) * nc1;    
    size_t nc3 = size_t(nc1) * nc1 * nc1;
    SearchIdxS sidx;
    for (size_t i = 0; i < nc3; i++) {
        int ix = i/nc2;
        int iy = (i%nc2)/nc1;
        int iz = i%nc1;
        ix = ix > nc1/2 ? ix - nc1 : ix;
        iy = iy > nc1/2 ? iy - nc1 : iy;
        iz = iz > nc1/2 ? iz - nc1 : iz;

        // compute r2_min        
        double dx = abs(ix) - 1 > 0 ? abs(ix) - 1 : 0;
        double dy = abs(iy) - 1 > 0 ? abs(iy) - 1 : 0;
        double dz = abs(iz) - 1 > 0 ? abs(iz) - 1 : 0;
        double r = sqrt(dx*dx + dy*dy + dz*dz) * cell_size;
        
        double min_rdi_sum = 2.0 * r / cutoff;
        if (min_rdi_sum > 2.0 * max_radius_) {
            continue;
        }

        // compute r2_max
        dx = abs(ix) + 1 > nc1 / 2 ? abs(ix) + 0.5 : abs(ix) + 1;
        dy = abs(iy) + 1 > nc1 / 2 ? abs(iy) + 0.5 : abs(iy) + 1;
        dz = abs(iz) + 1 > nc1 / 2 ? abs(iz) + 0.5 : abs(iz) + 1;
        r = sqrt(dx*dx + dy*dy + dz*dz) * cell_size;
        double max_rdi_sum = 2.0 * r / cutoff;
        
        // store search index
        sidx.dx = ix;
        sidx.dy = iy;
        sidx.dz = iz;
        sidx.max_rdi_sum = max_rdi_sum;
        sidx.min_rdi_sum = min_rdi_sum;

        search_idx_.push_back(sidx);
    }
}
    

int VerletListS::InitPairs(const double *rdi,
                           std::vector<std::vector<int> > *pairs)
{
    int npos = get_npos();
    double cutoff = get_cutoff();
    double max_radius = get_max_radius();
    double min_radius = get_min_radius();
    
    int nnz = 0;
    (*pairs).reserve(npos);
    for (int i = 0; i < npos; i++) {
        double t1 = (rdi[i] + max_radius)*cutoff/2.0 + max_radius;
        t1 = pow(t1, 3.0);
        double t2 = pow(rdi[i], 3.0);
        double t3 = pow(min_radius, 3.0);
        int size = (int)((t1 - t2)/t3 * 0.6) + 1;
        size = (size > npos ? npos : size);
        (*pairs)[i].reserve(size);
    }

    return nnz;
}


int VerletListS::FindPairs(const double *pos,
                           const double *rdi,
                           std::vector<std::vector<int> > *pairs)
{
    int nnz = 0;
    #pragma omp parallel
    {
        int len_sidx = search_idx_.size();
        int nc1 = get_nc1();
        int npos = get_npos();
        double cutoff = get_cutoff();
        double box_size = get_box_size();
        double max_radius = get_max_radius();
        double min_radius = get_min_radius();
        const std::vector<double> &pos0 = get_pos0();
        const std::vector<int> &head = get_head();
        const std::vector<int> &next = get_next();
        const std::vector<size_t> &cidx = get_cidx();
        size_t nc2 = nc1 * nc1;
        // for each cell list
        #pragma omp for reduction(+:nnz)
        for (int i = 0; i < npos; i++) {
            int idx1 = cidx[i];
            double radius_1 = rdi[i];
            double x1 = pos0[i * 3 + 0];
            double y1 = pos0[i * 3 + 1];
            double z1 = pos0[i * 3 + 2];
            int ix1 = idx1/nc2;
            int iy1 = (idx1%nc2)/nc1;
            int iz1 = idx1%nc1;
            (*pairs)[i].clear();
            // for each search index
            for (int k = 0; k < len_sidx; k++) {
                int dx = search_idx_[k].dx;
                int dy = search_idx_[k].dy;
                int dz = search_idx_[k].dz;
                double max_rdi_sum = search_idx_[k].max_rdi_sum;
                double min_rdi_sum = search_idx_[k].min_rdi_sum;
                int ix2 = ix1 + dx;
                int iy2 = iy1 + dy;
                int iz2 = iz1 + dz;
                ix2 = ix2 >= nc1 ? ix2 - nc1 : ix2;
                ix2 = ix2 < 0 ? ix2 + nc1 : ix2;
                iy2 = iy2 >= nc1 ? iy2 - nc1 : iy2;
                iy2 = iy2 < 0 ? iy2 + nc1 : iy2;
                iz2 = iz2 >= nc1 ? iz2 - nc1 : iz2;
                iz2 = iz2 < 0 ? iz2 + nc1 : iz2;
                size_t idx2 = ix2 * nc2 + iy2 * nc1 + iz2;                
                int j = head[idx2];
                if (radius_1 + min_radius > max_rdi_sum) {
                    while (j != -1) {
                        (*pairs)[i].push_back(j);
                        j = next[j];
                    }
                } else if (radius_1 + max_radius > min_rdi_sum) {
                    while (j != -1) {
                        double radius_2 = rdi[j];
                        if (radius_1 + radius_2 > max_rdi_sum) {
                            (*pairs)[i].push_back(j);
                        } else if (radius_1 + radius_2 > min_rdi_sum) {
                            double x2 = pos0[j * 3 + 0];
                            double y2 = pos0[j * 3 + 1];
                            double z2 = pos0[j * 3 + 2];
                            double rx = fabs(x2 - x1);
                            double ry = fabs(y2 - y1);
                            double rz = fabs(z2 - z1);
                            rx = (rx > box_size/2.0 ? rx - box_size : rx);
                            ry = (ry > box_size/2.0 ? ry - box_size : ry);
                            rz = (rz > box_size/2.0 ? rz - box_size : rz);
                            double rr = rx * rx + ry * ry + rz * rz;
                            if (rr == 0.0)
                                rr = 1.0e-10;
                            double s = 2.0 * sqrt(rr)/(radius_1 + radius_2);
                            if (s < cutoff) {
                                (*pairs)[i].push_back(j);
                            }
                        }
                        j = next[j];
                    }
                }
            }
            std::sort((*pairs)[i].begin(), (*pairs)[i].end());
            nnz += (*pairs)[i].size();
        }        
    }
    
    return nnz;
}

    
VerletListS::VerletListS(const int npos,
                         const double *rdi,
                         const double box_size,
                         const double cutoff)
    : VerletListBase(npos, rdi, box_size, cutoff)                     
{

}


VerletListS::~VerletListS()
{

}

} // namespace stokesdt