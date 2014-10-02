/**
 * @file   verlet_base.cc
 * @brief  VerletListBase class implementation
 */

#include <cmath>
#include <omp.h>
#include <algorithm>
#include "log.h"
#include "verlet_base.h"


namespace stokesdt {

VerletListBase::VerletListBase(const int npos,
                               const double *rdi,
                               const double box_size,
                               const double cutoff)
    : npos_(npos),
      box_size_(box_size),
      cutoff_(cutoff),
      rdi_(rdi, rdi + npos)
{

}


VerletListBase::~VerletListBase()
{

}


int VerletListBase::Init()
{
    if (npos_ <= 0) {
        LOG_ERROR("The specified number of paricles is less than"
                  "or equal to 0: %d.\n", npos_);
        return 0;
    }
    if (box_size_ <= 0.0){
        LOG_ERROR("The specified dimension of the simulation box is less than"
                  "or equal to 0.0: %g.\n", box_size_);
        return 0;        
    }
    if (cutoff_ <= 0.0){
        LOG_ERROR("The specified dimension of the simulation box is less than"
                  "or equal to 0.0: %g.\n", cutoff_);
        return 0;        
    }
    // get max and min radii
    max_radius_ = *std::max_element(rdi_.begin(), rdi_.end());
    min_radius_ = *std::min_element(rdi_.begin(), rdi_.end());
    if (min_radius_ <= 0.0) {
        LOG_ERROR("The specified particle radis is less than"
                  "or equal to 0.0: %g.\n", min_radius_);
        return 0;          
    }
    
    nc1_ = NumCells();    
    cell_size_ = box_size_/nc1_;
    size_t nc3 = size_t(nc1_) * nc1_ * nc1_;
    
    head_.reserve(nc3);
    next_.reserve(npos_);
    cidx_.reserve(npos_);
    pos0_.reserve(npos_ * 3);
    InitSearchIdx();

    int maxnnz = InitPairs(&rdi_[0], &pairs_);

    return maxnnz;
}


void VerletListBase::ConstructCellList(const double *pos,
                                       const double *rdi)
{
    #pragma omp parallel
    {        
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        size_t nc2 = (size_t)nc1_ * nc1_;
        size_t nc3 = (size_t)nc1_ * nc1_ * nc1_;
        size_t startidx = (nc3 + nthreads - 1)/nthreads * tid;
        size_t endidx = (nc3 + nthreads - 1)/nthreads * (tid + 1);
    
        // init head
        #pragma omp for
        for (size_t i = 0; i < nc3; i++) {
            head_[i] = -1;
        }

        // compute head
        for (int i = 0; i < npos_; i++) {
            double x = fmod(pos[3 * i + 0], box_size_);
            double y = fmod(pos[3 * i + 1], box_size_);
            double z = fmod(pos[3 * i + 2], box_size_);
            x = (x >= 0.0 ? x : box_size_ + x);
            y = (y >= 0.0 ? y : box_size_ + y);
            z = (z >= 0.0 ? z : box_size_ + z);
            int ix = (int)(x/box_size_ * nc1_);
            int iy = (int)(y/box_size_ * nc1_);
            int iz = (int)(z/box_size_ * nc1_);
            size_t idx = ix * nc2 + iy * nc1_ + iz;
            cidx_[i] = idx;
            pos0_[3 * i + 0] = x;
            pos0_[3 * i + 1] = y;
            pos0_[3 * i + 2] = z;            
            if (idx >= startidx && idx < endidx) {
                int h = head_[idx];
                head_[idx] = i;
                next_[i] = h;
            }
        }     
    }
}


int VerletListBase::Build(const double *pos)
{
    ConstructCellList(pos, &rdi_[0]);
    int nnz = FindPairs(pos, &rdi_[0], &pairs_);
    return nnz;
}


void VerletListBase::GetPairs(int *rowptr, int *colidx)
{
    rowptr[0] = 0;
    for (int i = 0; i < npos_; i++) {
        rowptr[i + 1] = rowptr[i] + pairs_[i].size();
    }
    #pragma omp parallel for
    for (int i = 0; i < npos_; i++) {
        int start = rowptr[i];
        int end = rowptr[i + 1];
        for (int j = start; j < end; j++) {
            colidx[j] = pairs_[i][j - start];
        }
    }
}

} // namespace stokesdt