/**
 * @file   mob_spme.cc
 * @brief  the MobSpme class implementation
 */

#include <vector>
#include <cmath>
#include <omp.h>
#include <mkl.h>
#include "ewald.h"
#include "mob_spme.h"
#include "log.h"
#include "matrix_io.h"


namespace stokesdt {

MobSpme::MobSpme(const int npos,
                 const double *rdi,
                 const double box_size,
                 const double xi,
                 const double rmax,
                 const int dim,
                 const int porder)
    : npos_(npos),
      box_size_(box_size),
      xi_(xi),
      rmax_(rmax),
      dim_(dim),
      porder_(porder),
      verlet_list_(npos, rdi, box_size, rmax)
{

}


MobSpme::~MobSpme()
{
    if (spme_ != NULL) {
        detail::DestroySpmeEngine(spme_);
    }
    if (real_mat_ != NULL) {
        detail::DestroySparseMatrix(real_mat_);
    }
}


bool MobSpme::Init()
{
    // check inputs
    if (npos_ <= 0) {
        LOG_ERROR("The specified number of particles is"
                  " less than or euqal to 0: %d.\n", npos_);
        return false;
    }
    if (box_size_ <= 0.0) {
        LOG_ERROR("The specified dimension of simulation box is"
                  " less than or euqal to 0.0: %lf.\n", box_size_);
        return false;
    }
    if (rmax_ <= 0.0) {
        LOG_ERROR("The specified real-space cutoff is less"
                  " than or euqal to 0.0: %lf.\n", rmax_);
        return false;
    }
    if (xi_ <= 0.0) {
        LOG_ERROR("The specified Ewald paramter is less than"
                  " or euqal to 0.0: %lf.\n", xi_);
        return false;
    }
    if (dim_ <= 0) {
        LOG_ERROR("The specified dimension of FFT mesh"
                  " less than or euqal to 0: %d.\n", dim_);
        return false;
    }    
    if (porder_ != 4 && porder_ != 6 && porder_ != 8) {
        LOG_ERROR("The specified interpolation order can only"
                  " be 4, 6 or 8: %d.\n", porder_);
        return false;
    }      
                 
    dim_mob_ = 3 * npos_;
    if (!detail::CreateSpmeEngine(npos_, box_size_, xi_,
                                  dim_, porder_, &spme_)) {
        LOG_ERROR("Creating Spme failed\n");
        return false;
    }

    int init_nnzb = verlet_list_.Init();

    if (!CreateSparseMatrix(npos_, init_nnzb, 3, &real_mat_)) {
        LOG_ERROR("Creating spare real matrix failed\n");
        return false;    
    }

    return true;
}


void MobSpme::BuildSparseReal(const double *pos, const double *rdi)
{    
    #pragma omp parallel
    {
        int *rowbptr = real_mat_->rowbptr;
        double *val = real_mat_->val;
        int *colbidx = real_mat_->colbidx;  
        double xi2 = xi_ * xi_;
        double xiaspi = xi_ / sqrt(M_PI);
        #pragma omp for
        for (int i = 0; i < npos_; i++) {
            double x1 = pos[3 * i + 0];
            double y1 = pos[3 * i + 1];
            double z1 = pos[3 * i + 2];
            double aa = rdi[i];
            double aa2 = aa * aa;
            double self_a = 1.0 / aa - xiaspi * (6.0 - 40.0 / 3.0 * xi2 * aa2);
            int start = rowbptr[i];
            int end = rowbptr[i + 1];
            for (int k = start; k < end; k++) {
                int j = colbidx[k];
                if (j != i) {
                    double x2 = pos[3 * j + 0];
                    double y2 = pos[3 * j + 1];
                    double z2 = pos[3 * j + 2];
                    double ab = rdi[j];
                    double ab2 = ab * ab;
                    double rx = drem(x2 - x1, box_size_);
                    double ry = drem(y2 - y1, box_size_);
                    double rz = drem(z2 - z1, box_size_);
                    double rr = rx * rx + ry * ry + rz * rz;
                    double r = sqrt(rr);
                    double xa;
                    double ya;
                    detail::scalars_ewald_F(xi_, r, aa, ab, &xa, &ya);
                    if (r < (aa + ab))  {
                        xa += -(1.5 - (aa2 + ab2) * 0.5 / rr) / r +
                                2.0 / (aa + ab) - 3.0 * r / 8.0 / (aa2 + ab2);
                        ya += -(0.75 + (aa2 + ab2) * 0.25 / rr) / r +
                                2.0 / (aa + ab) - 9.0 * r / 16.0 / (aa2 + ab2);
                    }
                    double ex = rx / r;
                    double ey = ry / r;
                    double ez = rz / r;
                    val[k * 9 + 0] = (xa - ya) * ex * ex + ya;
                    val[k * 9 + 1] = (xa - ya) * ex * ey;
                    val[k * 9 + 2] = (xa - ya) * ex * ez;
                    val[k * 9 + 3] = (xa - ya) * ex * ey;
                    val[k * 9 + 4] = (xa - ya) * ey * ey + ya;
                    val[k * 9 + 5] = (xa - ya) * ey * ez;
                    val[k * 9 + 6] = (xa - ya) * ex * ez;
                    val[k * 9 + 7] = (xa - ya) * ey * ez;
                    val[k * 9 + 8] = (xa - ya) * ez * ez + ya;
                } else {
                    val[k * 9 + 0] = self_a;
                    val[k * 9 + 1] = 0.0;
                    val[k * 9 + 2] = 0.0;
                    val[k * 9 + 3] = 0.0;
                    val[k * 9 + 4] = self_a;
                    val[k * 9 + 5] = 0.0;
                    val[k * 9 + 6] = 0.0;
                    val[k * 9 + 7] = 0.0;
                    val[k * 9 + 8] = self_a;
                }
            }
        } // #pragma omp for
    } // #pragma omp parallel
}


void MobSpme::Update(const double *pos, const double *rdi)
{
    detail::UpdateSpmeEngine(pos, rdi, spme_);
    int nnz = verlet_list_.Build(pos);
    ResizeSparseMatrix(nnz, real_mat_);
    verlet_list_.GetPairs(real_mat_->rowbptr, real_mat_->colbidx);
    BuildSparseReal(pos, rdi);
    int nm = dim();
    double *mat = (double *)malloc(sizeof(double) * nm * nm);
    SparseToDense(real_mat_, nm, mat);
}


void MobSpme::MulVector(const int num_rhs,
                        const double alpha,
                        const int ldvec_in,
                        const double *vec_in,
                        const double beta,
                        const int ldvec_out,
                        double *vec_out)
{
    detail::ComputeSpmeRecip(spme_, num_rhs,
                             alpha, ldvec_in, vec_in,
                             beta, ldvec_out, vec_out);
    SpMV3x3(real_mat_, num_rhs, alpha, ldvec_in, vec_in,
            1.0, ldvec_out, vec_out);
}


void MobSpme::RealMulVector(const int num_rhs,
                            const double alpha,
                            const int ldvec_in,
                            const double *vec_in,
                            const double beta,
                            const int ldvec_out,
                            double *vec_out)
{
    SpMV3x3(real_mat_, num_rhs, alpha, ldvec_in, vec_in,
            beta, ldvec_out, vec_out);
}


void MobSpme::RecipMulVector(const int num_rhs,
                             const double alpha,
                             const int ldvec_in,
                             const double *vec_in,
                             const double beta,
                             const int ldvec_out,
                             double *vec_out)
{
    detail::ComputeSpmeRecip(spme_, num_rhs,
                             alpha, ldvec_in, vec_in,
                             beta, ldvec_out, vec_out);
}

} // namespace stokesdt