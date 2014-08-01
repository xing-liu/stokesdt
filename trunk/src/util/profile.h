#ifndef PROFILE_H_
#define PROFILE_H_


#ifdef ENABLE_PROFILE_

#include <stdint.h>


namespace stokesdt {

namespace detail {

typedef enum {
    EWALD_TICKS              = 0,
    LUB_TICKS                = 1,
    CHOL_TICKS               = 2,
    CELLLIST_TICKS           = 3,
    REPULSION_TICKS          = 4,
    STOKES_NUM_TICKS
} StokesTicks_t;


const char ticks_name[STOKES_NUM_TICKS][128] = 
{
    "ewald",
    "lub",
    "cholesky",
    "celllist",
    "bond",
    "repulsion",
    "wrap",
    "interact",
    "interact2"
};


extern uint64_t stokes_etime[STOKES_NUM_TICKS];
extern uint64_t stokes_stime[STOKES_NUM_TICKS];


#define START_TIMER( idx )                                  \
        do {                                                \            
            stokes_stime[idx]  = __rdtsc();                 \
        } while ( 0 )
        
#define STOP_TIMER( idx )                                   \
        do {                                                \
            const uint64_t etime = __rdtsc();               \
            stokes_etime[idx] += etime - stokes_stime[idx]; \
        } while ( 0 )


inline void InitProfiler()
{
    uint64_t t = __rdtsc();
    sleep (1);
    stokes_freq = (double)(__rdtsc() - t);
    memset(stokes_etime, 0, sizeof(uint64_t) * STOKES_NUM_TICKS);
    LOG(3, "CPU Freq = %.3lf GHz", stokes_freq/1e9);
}


inline void ResetProfiler()
{
    memset(stokes_etime, 0, sizeof(uint64_t) * STOKES_NUM_TICKS);
}


inline void PrintProfiler()
{
    for (int k = 0; k < STOKES_NUM_TICKS; k++) {
        if (stokes_etime[k] != 0) {
            LOG(2, "%12s", ticks_name[k]);
            LOG(2, ",\t%.3le secs", (double)stokes_etime[k]/stokes_freq);
            LOG(2, "\n");
        }
    }
}

} // namespace detail
 
} // namespace stokesdt

#endif // #ifdef ENABLE_PROFILE_


#endif //  PROFILE_H_
