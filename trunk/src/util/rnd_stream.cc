/**
 *  @file   rnd_stream.cc
 *  @brief  RndStream class implementation
 */

#include <mkl_vsl.h>
#include "rnd_stream.h"
#include "log.h"
#include "profile.h"


namespace stokesdt {

RndStream::RndStream(const unsigned int seed) : seed_(seed)
{

}


bool RndStream::Init()
{
    LOG(3, "\n        Initializes RndStream\n");
    LOG(3, "        ---------------------\n");
    LOG(3, "Sedd = %d\n", seed_);

    VSLStreamStatePtr vslstream;
    if (vslNewStream(&vslstream, VSL_BRNG_SFMT19937, seed_) != VSL_STATUS_OK) {
        return false;
    } else {
        stream_ = vslstream;
    }

    return true;
}


RndStream::~RndStream()
{
    VSLStreamStatePtr vslstream = (VSLStreamStatePtr)stream_;
    vslDeleteStream(&vslstream);
}


void RndStream::Gaussian(const double mean,
                         const double sigma,
                         const int num_rhs,
                         const int len_v,
                         const int ldv,
                         double *v)
{
    START_TIMER(detail::RANDOM_TICKS);
    
    VSLStreamStatePtr vslstream = (VSLStreamStatePtr)stream_;
    if (ldv == len_v) {
        vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                      vslstream, num_rhs * len_v, v, mean, sigma);
    } else {
        for (int i = 0; i < num_rhs; i++) {
            vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                          vslstream, len_v, &(v[i * ldv]), mean, sigma);
        }
    }

    STOP_TIMER(detail::RANDOM_TICKS);
}


void RndStream::Uniform(const double low,
                        const double high,
                        const int num_rhs,
                        const int len_v,
                        const int ldv,
                        double *v)
{
    VSLStreamStatePtr vslstream = (VSLStreamStatePtr)stream_;
    if (ldv == len_v) {
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE,
                     vslstream, num_rhs * len_v, v, low, high);
    } else {
        for (int i = 0; i < num_rhs; i++) {
            vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE,
                         vslstream, len_v, &(v[i * ldv]), low, high);
        }
    }
}

} // namespace stokesdt
