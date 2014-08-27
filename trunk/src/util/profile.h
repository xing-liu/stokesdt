/**
 * @file   profile.h
 * @brief  Profiler definition
 */

#ifndef PROFILE_H_
#define PROFILE_H_


#include <stdint.h>


namespace stokesdt {

namespace detail {

typedef enum {
    MOB_TICKS                = 0,
    BRWN_TICKS               = 1,
    STERIC_TICKS             = 2,
    BONDED_TICKS             = 3,
    RANDOM_TICKS             = 4,
    PARTICLE_TICKS           = 5,
    STOKES_NUM_TICKS
} StokesTicks_t;


const char ticks_name[STOKES_NUM_TICKS][128] = 
{
    "Mob",
    "Brownian",
    "Steric",
    "Bonded",
    "Random",
    "Particle"
};

extern uint64_t stokes_etime[STOKES_NUM_TICKS];
extern uint64_t stokes_stime[STOKES_NUM_TICKS];

void InitProfiler();

void ResetProfiler();

void PrintProfiler();

} // namespace detail
 
} // namespace stokesdt


#ifdef ENABLE_PROFILE_

#define START_TIMER(idx)                                      \
        do {                                                  \            
            stokesdt::detail::stokes_stime[idx]  = __rdtsc(); \
        } while ( 0 )
        
#define STOP_TIMER(idx)                                      \
        do {                                                 \
            const uint64_t etime = __rdtsc();                \
            stokesdt::detail::stokes_etime[idx] +=           \
                etime - stokesdt::detail::stokes_stime[idx]; \
        } while ( 0 )

#else

#define START_TIMER(idx)
        
#define STOP_TIMER(idx)
        
#endif // #ifdef ENABLE_PROFILE_


#endif //  PROFILE_H_
