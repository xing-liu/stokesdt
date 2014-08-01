/**
 * @file   force_base.cc
 * @brief  ForceBase class implementation
 */

#include <string.h>
#include "force_base.h"


namespace stokesdt {

ForceBase::ForceBase(const int npos) : npos_(npos)
{

}


ForceBase::~ForceBase()
{

}


void ForceBase::Compute(const double *pos, const double *rdi, double *f)
{
    memset(f, 0, sizeof(double) * npos_ * 3);
    Accumulate(pos, rdi, f);
}

} // namespace stokesdt