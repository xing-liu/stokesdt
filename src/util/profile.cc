/**
 * @file   profile.cc
 * @brief  Profiler implementation
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "profile.h"
#include "log.h"


namespace stokesdt {

namespace detail {
    
uint64_t stokes_etime[STOKES_NUM_TICKS];
uint64_t stokes_stime[STOKES_NUM_TICKS];
static double stokes_freq;

void InitProfiler()
{
    uint64_t t = __rdtsc();
    sleep (1);
    stokes_freq = (double)(__rdtsc() - t);
    memset(stokes_etime, 0, sizeof(uint64_t) * STOKES_NUM_TICKS);
}


void ResetProfiler()
{
    memset(stokes_etime, 0, sizeof(uint64_t) * STOKES_NUM_TICKS);
}


void PrintProfiler()
{
    LOG(3, "\n        Profiling Results\n");
    LOG(3, "        -----------------\n");
    LOG(3, "        CPU Freq = %.3lf GHz\n\n", stokes_freq/1e9);

    double total_time = 0.0;
    for (int k = 0; k < STOKES_NUM_TICKS; k++) {
        if (stokes_etime[k] != 0.0) {
            LOG(3, "%12s", ticks_name[k]);
            LOG(3, ",\t%.3g secs", (double)stokes_etime[k]/stokes_freq);
            LOG(3, "\n");
            total_time += (double)stokes_etime[k]/stokes_freq;
        }
    }
    LOG(3, "%12s,\t%.3g secs\n", "Total", total_time);
}

} // namespace detail
 
} // namespace stokesdt
