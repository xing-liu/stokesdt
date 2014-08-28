#include <cstdio>
#include <getopt.h>
#include <cstring>
#include <sys/time.h>

#include "stokes.h"


using namespace stokesdt;


const struct option long_options[] = {
    {"help",    0, NULL, 'h'},
    {"version", 0, NULL, 'v'},
    {"nsteps",  1, NULL, 'n'},
    {"config",  1, NULL, 10},
    {"model",   1, NULL, 20},
    {"xyz",     1, NULL, 30},
    {"log",     1, NULL, 40},
    {NULL,      0, NULL,  0}
};

const char *const short_options = ":hvn:";
const char *version_info = "0.1.0";


static void Usage(char *call)
{
    fprintf(stderr, "Usage: %s [OPTIONS]\n", call);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-h or --help         Display this information\n");
    fprintf(stderr, "\t-v or --version      Display version information\n");
    fprintf(stderr, "\t--config             Config file\n");
    fprintf(stderr, "\t--model              Model file\n");
    fprintf(stderr, "\t--xyz                XYZ file\n");
    fprintf(stderr, "\t--log                log file\n");
    fprintf(stderr, "\t-n or --nsteps       Number of steps\n");
}


static void PrintVersion(char *call)
{
    fprintf(stdout, "%s version %s\n", call, version_info);
}


int main(int argc, char **argv)
{
    char *config_file = NULL;
    char *model_file = NULL;
    char *xyz_file = NULL;
    char *log_str = NULL;
    int num_steps = 1;
    
    // parse arguments
    int c = 0;
    while ((c = getopt_long(argc, argv, short_options,
                            long_options, NULL)) != -1) {
        switch (c) {
            case 'h':
                Usage(argv[0]);
                return 0;
            case 'v':
                PrintVersion(argv[0]);
                return 0;
            case 'n':
                num_steps = atoi(optarg);
                if (num_steps <= 0) {
                    fprintf(stderr, "Invalid argument --nsteps %d\n", num_steps);
                }
                break;
                return 0;
            case 10:
                config_file = strdup(optarg);
                break;           
            case 20:
                model_file = strdup(optarg);
                break;           
            case 30:
                xyz_file = strdup(optarg);
                break;
            case 40:
                log_str = strdup(optarg);
                break;                
            case ':':
                fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                return -1;
            case '?':
                fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                return -1;
            default:
                Usage(argv[0]);
                return -1;
        }
    }
    if (NULL == model_file ||
        NULL == xyz_file   ||
        NULL == config_file) {
        Usage(argv[0]);
        return -1;
    }

    const int kLenBlk = 100;
    struct timeval tv1;
    struct timeval tv2;

    // Start logging
    if (NULL != log_str) {
        char *log_file = strtok(log_str, ":"); 
        if (log_file != NULL &&
            !InitLogger(log_file)) {
            fprintf(stderr, "Failed to start logging: %s\n", log_file);
            return false;
        } else {
            char *level_str = strtok(NULL, ":");
            int log_level = 1;
            if (level_str != NULL) {
                log_level = atoi(level_str);
            } else {
                log_level = 3;
            }
            if (log_level >=0) {
                SetLoggerLevel(log_level);
            }
        }
    }
    
    // Create simulation
    BDSimBox *bd_sim = new BDSimBox(config_file, model_file, xyz_file);
    if (!bd_sim->Init()) {
        fprintf(stderr, "Failed to initialize BDSimBox\n");
        return -1;
    }

    fprintf(stdout, "Start BD simulation\n");
    fprintf(stdout, "-------------------\n");
    
    // Start simulation
    double totaltime = 0.0;   
    for (int i = 0; i < num_steps; i+=kLenBlk) {
        int nsteps_blk  = (kLenBlk + i >= num_steps ? num_steps - i : kLenBlk);
        gettimeofday(&tv1, NULL);
        
        bd_sim->Advance(nsteps_blk);

        gettimeofday(&tv2, NULL);
        timersub(&tv2, &tv1, &tv1);
        totaltime += tv1.tv_sec + tv1.tv_usec/1e6;
        
        fprintf(stdout, "finished steps: %10d (%7.3f%%),"
                        "    elapsed time: %10.2f secs\r",
                        i + nsteps_blk,
                        100.0 * (i + nsteps_blk)/num_steps,
                        totaltime);
        fflush(stdout);
    }
    fprintf(stdout, "\n\nBD simulation finished at %.3g secs, %.3g secs/step\n",
        totaltime, totaltime/num_steps);
       
    // clean up
    delete bd_sim;
    free(model_file);
    free(xyz_file);
    free(config_file);
    free(log_str);
    DestroyLogger();
    
    return 0;
}
