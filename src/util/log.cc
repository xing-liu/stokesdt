/**
 * @file   log.cc
 * @brief  Logger implementation
 */
 
#ifdef ENABLE_LOG_

#include "log.h"


/// Log level
int log_level = 3;
/// Log file
FILE *fp_log = NULL;


/** @brief Initialize the logger
 *
 *  @param[in] file  log file name 
 */
bool InitLogger(char *file)
{
    fp_log = fopen(file, "w+");
    if (NULL == fp_log) {
        fprintf(stderr, "Can't open %s\n", file);
        return false;
    }

    return true;
}


/** @brief Destroy the logger
 */
void DeinitLogger()
{
    if (fp_log != NULL) {
        fclose(fp_log);
    }
}

#endif // ENABLE_LOG_ 
