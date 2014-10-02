/**
 * @file   force_steric.cc
 * @brief  StericForce class implementation
 */

#include <cmath>
#include "force_steric.h"
#include "log.h"
#include "profile.h"


namespace stokesdt {

void StericForce::ForceKernel(const double *pos,
                              const double *rdi, 
                              double *f)
{
    int npos = get_npos();
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < npos; i++)
        {
            double x1 = pos[3 * i + 0];
            double y1 = pos[3 * i + 1];
            double z1 = pos[3 * i + 2];
            double aa = rdi[i];
            int start = rowptr_[i];
            int end = rowptr_[i + 1];
            for (int k = start; k < end; k++)
            {
                int j = colidx_[k];
                if (i == j)
                    continue;
                double x2 = pos[3 * j + 0];
                double y2 = pos[3 * j + 1];
                double z2 = pos[3 * j + 2];
                double ab = rdi[j];
                double rx = drem(x1 - x2, box_size_);
                double ry = drem(y1 - y2, box_size_);
                double rz = drem(z1 - z2, box_size_);
                double rr = rx * rx + ry * ry + rz * rz;
                double r = sqrt(rr);
                double ex = rx / r;
                double ey = ry / r;
                double ez = rz / r;
                double force =
                    steric_k0_ * ((aa + ab) * steric_r0_ / 2.0 - r);
                f[3 * i + 0] +=  force * ex;
                f[3 * i + 1] +=  force * ey;
                f[3 * i + 2] +=  force * ez;
            }
        }
    }
}


StericForce::StericForce(const int npos,
                         const double *rdi,
                         const double box_size, 
                         const double steric_r0,
                         const double steric_k0)
    : ForceBase(npos),
      steric_k0_(steric_k0),
      steric_r0_(steric_r0),
      box_size_(box_size),
      rowptr_(npos + 1),
      verlet_list_(npos, rdi, box_size, steric_r0_)
{

}


StericForce::~StericForce()
{
}


bool StericForce::Init()
{
    if (steric_k0_ <= 0.0) {
        LOG_ERROR("The steric coefficient is less than or equal to 0.0\n");
        return false;
    }
    LOG(3, "\n        Initializes StericForce\n");
    LOG(3, "        -----------------------\n");    
    LOG(3, "Steric-r0 = %g\n", steric_r0_);
    LOG(3, "Steric-k0 = %g\n", steric_k0_);
    
    int nnz = verlet_list_.Init();
    colidx_.reserve(nnz);

    return true;
}


void StericForce::Accumulate(const double *pos, const double *rdi, double *f)
{
   START_TIMER(detail::STERIC_TICKS);
   
   // update Verlet list
   int nnz = verlet_list_.Build(pos);
   colidx_.reserve(nnz);
   verlet_list_.GetPairs(&rowptr_[0], &colidx_[0]);

   // comput steric forces
   ForceKernel(pos, rdi, f);

   STOP_TIMER(detail::STERIC_TICKS);
}

} // namespace stokesdt
