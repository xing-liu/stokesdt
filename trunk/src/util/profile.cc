#ifdef ENABLE_PROFILE_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "profile.h"


namespace stokesdt {

namespace detail {
    
uint64_t stokes_etime[STOKES_NUM_TICKS];
uint64_t stokes_stime[STOKES_NUM_TICKS];
static double stokes_freq;

} // namespace detail
 
} // namespace stokesdt


#endif // #ifdef ENABLE_PROFILE_
