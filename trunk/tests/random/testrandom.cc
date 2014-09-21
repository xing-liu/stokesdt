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
    {"nrhs",    1, NULL, 40},
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
    fprintf(stderr, "\t--nrhs               Number of right-hand sides\n");
}

int main(int argc, char **argv)
{
    char *model_file = NULL;
    char *xyz_file = NULL;
    int num_vecs = 5;
    int num_rhs = 10;
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
            case 40:
                num_rhs = atoi(optarg);
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
    int blocksize = 1000;
    double *z = (double *)malloc(sizeof(double) * ldm * blocksize);
    assert(z != NULL);
    RndStream *rnd = new RndStream(1234);
    assert(rnd->Init());
    double *y = (double *)malloc(sizeof(double) * ldm * blocksize);
    assert(y != NULL);
    double *cov = (double *)malloc(sizeof(double) * nm * ldm);
    assert(cov != NULL);
    memset(cov, 0, sizeof(double) * nm * ldm);
    // create randomized
    BrwnRandom *random = new BrwnRandom(nm, num_vecs, blocksize);
    assert(random->Init());
    // init
    random->Compute(mob, 1, ldm, z, ldm, y);

    // compute randomized sampling
    MatrixIO *matrix_io = new MatrixIO();
    for (int i = 0; i < num_rhs; i+=blocksize) {
        int kend = i + blocksize <= num_rhs ? blocksize : num_rhs - i;
        // compute vectors
        rnd->Gaussian(0.0, 1.0, kend, nm, ldm, z);
        random->Compute(NULL, kend, ldm, z, ldm, y);
        // statistics
        for (int k = 0; k < kend; k++) {
            // update cov
            #pragma omp parallel for
            for (int p = 0; p < nm; p++) {
                #pragma simd
                for (int q = 0; q < nm; q++) {
                    cov[p * ldm + q] += y[k * ldm + p] * y[k * ldm + q];
                }
            }         
        }
        #if 0
        double scalar = (double)i + kend;
        printf("i = %d %d %g\n", i, kend, scalar);
        double error;
        if (scalar != 0.0) {
            cblas_dscal(nm * ldm, 1.0/scalar, cov, 1);
        }
        matrix_io->Compare(nm, nm, ldm, mat, ldm, cov, &error);
        if (scalar != 0) {
            cblas_dscal(nm * ldm, scalar, cov, 1);
        }
        printf("error = %g\n", error);
        #endif
    }
    double scalar = (double)num_rhs;
    double error;
    if (scalar != 0.0) {
        cblas_dscal(nm * ldm, 1.0/scalar, cov, 1);
    }
    matrix_io->Compare(nm, nm, ldm, mat, ldm, cov, &error);
    printf("error = %g\n", error);
    
    // clean up
    delete parser;
    delete mob;
    delete random;
    delete rnd;
    delete matrix_io;
    free(cov);
    free(y);
    free(z);
    free(pos);
    free(rdi);
    free(mat);
    free(model_file);
    free(xyz_file);
    
    return 0;
}
