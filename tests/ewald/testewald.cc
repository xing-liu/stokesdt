#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>

#include "mob_ewald.h"
#include "molecule_io.h"
#include "matrix_io.h"


using namespace stokesdt;

const struct option long_options[] = {
    {"model",   1, NULL, 10},
    {"xyz",     1, NULL, 20},
    {"ref",     1, NULL, 30},
    {NULL, 0, NULL, 0}
};

const char *const short_options = ":h";

static void usage(char *call)
{
    fprintf(stderr, "Usage: %s [OPTIONS]\n", call);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-h or --help         Display this information\n");
    fprintf(stderr, "\t--model              Model file\n");
    fprintf(stderr, "\t--xyz                XYZ file\n");
    fprintf(stderr, "\t--ref                Reference file\n");
}
int main(int argc, char **argv)
{
    char *model_file = NULL;
    char *xyz_file = NULL;
    char *ref_file = NULL;
    /* parse arguments */
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
                ref_file = strdup(optarg);
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
        NULL == xyz_file ||
        NULL == ref_file) {
        usage(argv[0]);
        return -1;
    }
    fprintf(stdout, "\nTest Ewald summation\n");
    fprintf(stdout, "--------------------\n");
    
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

    // compute mobility
    MobEwald *mob = new MobEwald(npos, rdi, box_size, 1.0e-12);
    assert(mob->Init());
    int nm = mob->dim();
    int ldm = nm;
    mob->Update(pos, rdi);    
    double *mat = (double *)malloc(sizeof(double) * nm * ldm);
    assert(mat != NULL);    
    mob->GetMatrix(ldm, mat);

    // check results
    MatrixIO *matrix_io = new MatrixIO();
    double *mat0 = (double *)malloc(sizeof(double) * nm * ldm);
    assert(mat0 != NULL);
    matrix_io->Read(ref_file, nm, nm, ldm, mat0);
    matrix_io->Write(nm, nm, ldm, mat, "mob.dat");
    double error;
    matrix_io->Compare(nm, nm, ldm, mat, ldm, mat0, &error);
    printf("error = %le\n", error);
    if (error < 1.0e-6) {
        printf("PASS.\n");
    } else {
        printf("FAILED.\n\n");
    }
    
    // clean up
    delete matrix_io;
    delete parser;
    delete mob;
    free(pos);
    free(rdi);
    free(mat);
    free(mat0);
    free(model_file);
    free(ref_file);
    free(xyz_file);
    
    return 0;
}
