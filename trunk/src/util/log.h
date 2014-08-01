/**
 * @file   log.h
 * @brief  Defines the logger
 */

#ifndef LOG_H_
#define LOG_H_


#ifdef ENABLE_LOG_


#include <string>

class Logger{
  public:
    static Logger* Instance();
    bool openLogFile(std::string logFile);
    void writeToLogFile();
    bool closeLogFile();
 private:
   Logger(){};  // Private so that it can  not be called
   Logger(Logger const&){};             // copy constructor is private
   Logger& operator=(Logger const&){};  // assignment operator is private
   static Logger* m_pInstance;
};


#include <stddef.h>
#include "logger.h"

// Global static pointer used to ensure a single instance of the class.
Logger* Logger::m_pInstance = NULL; 
 
Logger* Logger::Instance()
{
   if (!m_pInstance)   // Only allow one instance of class to be generated.
      m_pInstance = new Logger;
 
   return m_pInstance;
}

bool Logger::openLogFile(std::string _logFile)
{

}


#include <cstdio>

extern int log_level;
extern FILE *fp_log;

/// Log a message
#define LOG(level, fmt, args...)                             \
        do {                                                 \
            if (level <= log_level && fp_log != NULL) {      \
                fprintf(fp_log, fmt, ##args);                \
                fflush(fp_log);                              \
            }                                                \
        } while (0)
        
/// Log an error        
#define LOG_ERROR(fmt, args...)                              \
        do {                                                 \
            if (fp_log != NULL) {                            \
                fprintf(fp_log, "%s:%d: %s(): "fmt, __FILE__,\
                    __LINE__, __func__, ##args);             \
                fflush(fp_log);                              \
            }                                                \
            fprintf(stderr, "%s:%d: %s(): "fmt, __FILE__,    \
                __LINE__, __func__, ##args);                 \
            fflush(stderr);                                  \
        } while (0)

bool InitLogger(char *file);

void DeinitLogger();            
        
#else
/// Log a message
#define LOG(level, fmt, args...)

/// Log an error 
#define LOG_ERROR(fmt, args...)                              \
        do {                                                 \
            fprintf(stderr, "%s:%d: %s(): "fmt, __FILE__,    \
                    __LINE__, __func__, ##args);             \
            fflush(stderr);                                  \
        } while (0)
        
#endif // #ifdef ENABLE_LOG_


#define LOG_WARN(fmt, args...)  LOG(1, fmt, ##args); 

        
#endif // LOG_H_
