#include <cstdio>
#include <getopt.h>
#include <cstring>

#include "stokes.h"


using namespace stokesdt;


const struct option long_options[] = {
    {"config",  1, NULL, 10},
    {"model",   1, NULL, 20},
    {"xyz",     1, NULL, 30},
    {"nsteps",  1, NULL, 40},
    {NULL,      0, NULL,  0}
};

const char *const short_options = ":h";

static void usage (char *call)
{
    fprintf (stderr, "Usage: %s [OPTIONS]\n", call);
    fprintf (stderr, "Options:\n");
    fprintf (stderr, "\t-h or --help         Display this information\n");
    fprintf (stderr, "\t--config             config file\n");
    fprintf (stderr, "\t--model              Model file\n");
    fprintf (stderr, "\t--xyz                XYZ file\n");
    fprintf (stderr, "\t--nsteps             Number of steps\n");
}


int main(int argc, char **argv)
{
    char *config_file = NULL;
    char *model_file = NULL;
    char *xyz_file = NULL;
    int nsteps = 1;
    
    // parse arguments
    int c = 0;
    while ((c = getopt_long(argc, argv, short_options,
                            long_options, NULL)) != -1) {
        switch (c) {
            case 'h':
                usage(argv[0]);
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
                nsteps = atoi(optarg);
                if (nsteps <= 0) {
                    fprintf(stderr, "Invalid argument --nsteps %d\n", nsteps);
                }
                break;
            case ':':
                fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                return -1;
            case '?':
                fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                return -1;
            default:
                usage (argv[0]);
                return -1;
        }
    }
    if (NULL == model_file ||
        NULL == xyz_file   ||
        NULL == config_file) {
        usage(argv[0]);
        return -1;
    }
    
    // Start simulation
    BDSimBox *bd_sim = new BDSimBox(config_file, model_file, xyz_file);
    if (!bd_sim->Init()) {
        fprintf(stderr, "Initialize BDSimBox failed\n");
        return -1;
    }

    bd_sim->Advance(nsteps);
    
    
    // clean up
    delete bd_sim;
    free(model_file);
    free(xyz_file);
    free(config_file);
    
    return 0;
}
