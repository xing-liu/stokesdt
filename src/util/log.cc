/**
 * @file   log.cc
 * @brief  Logger implementation
 */
 

#include "log.h"


namespace stokesdt {

namespace detail {

/// the log level
int log_level = 3;
/// the log file
FILE *fp_log = NULL;

} // namespace detail


bool InitLogger(char *file)
{
    detail::fp_log = fopen(file, "w+");
    if (NULL == detail::fp_log) {
        fprintf(stderr, "Can't open %s\n", file);
        return false;
    }

    return true;
}


void DestroyLogger()
{
    if (detail::fp_log != NULL) {
        fclose(detail::fp_log);
    }
}


void SetLoggerLevel(int level)
{
    detail::log_level = level;
}


} // namespace stokesdt
