#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <mkl.h>
#include <string.h>

#include "stokes_mob.h"
#include "stokes_util.h"
#include "stokes_brwn.h"


using namespace stokesdt;

const struct option long_options[] = {
    {"model",   1, NULL, 10},
    {"xyz",     1, NULL, 20},
    {"nvecs",   1, NULL, 30},
    {NULL,      0, NULL,  0}
};

const char *const short_options = ":h";

static void usage (char *call)
{
    fprintf(stderr, "Usage: %s [OPTIONS]\n", call);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-h or --help         Display this information\n");
    fprintf(stderr, "\t--model              Model file\n");
    fprintf(stderr, "\t--xyz                XYZ file\n");
    fprintf(stderr, "\t--nvecs              Number of random vectors\n");
}

int main(int argc, char **argv)
{
    char *model_file = NULL;
    char *xyz_file = NULL;
    int num_vecs = 5;
    // parse arguments
    int c = 0;
    while ((c = getopt_long(argc, argv, short_options,
                            long_options, NULL)) != -1) {
        switch (c) {
            case 'h':
                usage(argv[0]);
                return 0;
            case 10:
                model_file = strdup(optarg);
                break;           
            case 20:
                xyz_file = strdup(optarg);
                break;           
            case 30:
                num_vecs = atoi(optarg);
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
        NULL == xyz_file) {
        usage(argv[0]);
        return -1;
    }
    
    fprintf(stdout, "\nTest Random\n");
    fprintf(stdout, "-------------\n");
    
    // input particles
    MoleculeIO * parser = new MoleculeIO();
    assert(parser->ParseModel(model_file, " \t", "#"));
    int npos;
    double box_size;
    assert(parser->ParseXYZ(xyz_file, &npos, &box_size));
    double *pos = (double *)malloc(sizeof(double) * npos * 3);
    double *rdi = (double *)malloc(sizeof(double) * npos);
    assert(pos != NULL && rdi != NULL);
    parser->GetParticles(pos, rdi);

    // compute mobility matrix
    MobEwald *mob = new MobEwald(npos, rdi, box_size, 1.0e-12);
    assert(mob->Init());
    int nm = mob->dim();
    int ldm = nm;
    mob->Update(pos, rdi);    
    double *mat = (double *)malloc(sizeof(double) * nm * ldm);
    assert(mat != NULL);    
    mob->GetMatrix(ldm, mat);

    // create buffers
    int num_rhs = 3;
    double *z = (double *)malloc(sizeof(double) * ldm * num_rhs);
    assert(z != NULL);
    RndStream *rnd = new RndStream(1234);
    assert(rnd->Init());
    rnd->Gaussian(0.0, 1.0, num_rhs, nm, ldm, z);
    double *y = (double *)malloc(sizeof(double) * ldm * num_rhs);
    assert(y != NULL);

    // compute randomized
    BrwnRandom *random = new BrwnRandom(nm, num_vecs, num_rhs);
    assert(random->Init());

    random->Compute(mob, num_rhs, ldm, z, ldm, y);
    
    // clean up
    delete parser;
    delete mob;
    delete random;
    delete rnd;
    free(y);
    free(z);
    free(pos);
    free(rdi);
    free(mat);
    free(model_file);
    free(xyz_file);
    
    return 0;
}
