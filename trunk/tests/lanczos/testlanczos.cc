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
}

int main(int argc, char **argv)
{
    char *model_file = NULL;
    char *xyz_file = NULL;   
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
    
    fprintf(stdout, "\nTest Lanczos\n");
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
    rnd->Uniform(0.0, 1.0, num_rhs, nm, ldm, z);
    double *y = (double *)malloc(sizeof(double) * ldm * num_rhs);
    double *y0 = (double *)malloc(sizeof(double) * ldm * num_rhs);
    assert(y != NULL && y0 != NULL);

    // compute lanczos
    int maxiters = 100;
    int maxnrhs = 16;
    BrwnLanczos *lanczos = new BrwnLanczos(nm, maxiters, maxnrhs, 1.0e-7);
    assert(lanczos->Init());
    
    // eigenvalue, mat2 = mat^1/2
    double *v = (double *)malloc(sizeof(double) * nm * ldm);
    assert(v != NULL);
    double *mat2 = (double *)malloc(sizeof(double) * nm * ldm);
    assert(mat2 != NULL);
    double *w = (double *)malloc(sizeof(double) * nm);
    assert(w != NULL);
    MatrixIO *matrix_io = new MatrixIO();
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nm, mat, ldm, w);
    memcpy(v, mat, sizeof(double) * nm * ldm);
    for(int i = 0; i < nm; i++) {
        cblas_dscal(nm, sqrt(w[i]), &mat[i], ldm);
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                nm, nm, nm, 1.0, mat, ldm, v, ldm,
                0.0, mat2, ldm);

    // check results
    num_rhs = 1;
    lanczos->Compute(mob, num_rhs, ldm, z, ldm, y);    
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                nm, num_rhs, nm, 1.0, mat2, ldm, z, ldm,
                0.0, y0, ldm);
    double error1;
    matrix_io->Compare(num_rhs, nm, ldm, y, ldm, y0, &error1);
    printf("error 1 = %le\n", error1);
    
    num_rhs = 3;
    lanczos->Compute(mob, num_rhs, ldm, z, ldm, y);    
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                nm, num_rhs, nm, 1.0, mat2, ldm, z, ldm,
                0.0, y0, ldm);
    double error3;
    matrix_io->Compare(num_rhs, nm, ldm, y, ldm, y0, &error3);
    printf("error 3 = %le\n", error3);
        

    if (error1 < 1.0e-6 && error3 < 1.0e-6) {
        printf("PASS.\n");
    } else {
        printf("FAILED.\n\n");
    }
    
    // clean up
    delete matrix_io;
    delete parser;
    delete mob;
    delete lanczos;
    delete rnd;
    free(y);
    free(y0);
    free(mat2);
    free(v);
    free(z);
    free(pos);
    free(rdi);
    free(mat);
    free(w);
    free(model_file);
    free(xyz_file);
    
    return 0;
}
