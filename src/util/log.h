/**
 * @file   log.h
 * @brief  Defines the logger
 */

#ifndef LOG_H_
#define LOG_H_


#include <cstdio>


namespace stokesdt {

namespace detail {

/// the log level
extern int log_level;
/// the log file
extern FILE *fp_log;

} // namespace detail

/** 
 * @brief   Initializes the logger
 * @param[in] file  the log file name 
 */
bool InitLogger(char *file);

/// @brief  Destroys the logger
void DestroyLogger();

/// @brief  Sets the log level
void SetLoggerLevel(int level);


} // namespace stokesdt


/// Write a message into the log file
#define LOG(level, fmt, args...)                                  \
        do {                                                      \
            if (level <= stokesdt::detail::log_level &&           \
                stokesdt::detail::fp_log != NULL) {               \
                fprintf(stokesdt::detail::fp_log, fmt, ##args);   \
                fflush(stokesdt::detail::fp_log);                 \
            }                                                     \
        } while (0)
        
/// Write an error message into the log file      
#define LOG_ERROR(fmt, args...)                                   \
        do {                                                      \
            if (stokesdt::detail::fp_log != NULL) {               \
                fprintf(stokesdt::detail::fp_log,                 \
                        "ERROR: %s:%d: %s(): "fmt, __FILE__,      \
                        __LINE__, __func__, ##args);              \
                fflush(stokesdt::detail::fp_log);                 \
            }                                                     \
            fprintf(stderr, "ERROR: %s:%d: %s(): "fmt, __FILE__,  \
                __LINE__, __func__, ##args);                      \
            fflush(stderr);                                       \
        } while (0)

/// Write an warning message into the log file   
#define LOG_WARN(fmt, args...)  LOG(1, "WARNING: "fmt, ##args);

        
#endif // LOG_H_
