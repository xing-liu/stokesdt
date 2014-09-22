#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <cstring>

#include "stokes.h"
#include "mob_debug.h"


using namespace stokesdt;

const struct option long_options[] = {
    {"model",   1, NULL, 10},
    {"xyz",     1, NULL, 20},
    {"xi",      1, NULL, 30},
    {"rmax",    1, NULL, 40},
    {"dim",     1, NULL, 50},
    {"porder",  1, NULL, 60},
    {NULL, 0, NULL, 0}
};

const char *const short_options = ":h";

static void usage (char *call)
{
    fprintf(stderr, "Usage: %s [OPTIONS]\n", call);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-h or --help         Display this information\n");
    fprintf(stderr, "\t--model              Model file\n");
    fprintf(stderr, "\t--xyz                XYZ file\n");
    fprintf(stderr, "\t--xi                 Ewald paramter\n");    
    fprintf(stderr, "\t--rmax               Real-space cutoff\n");
    fprintf(stderr, "\t--dim                Dimension of FFT grid\n");
    fprintf(stderr, "\t--porder             Interpolation order\n");  
}


int main(int argc, char **argv)
{
    char *model_file = NULL;
    char *xyz_file = NULL;
    double xi = 0.5;
    double rmax = 4.0;
    int dim = 64;
    int porder = 4;
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
                xi = atof(optarg);
                break;                
            case 40:
                rmax = atof(optarg);
                break;                
            case 50:
                dim = atoi(optarg);
                break;                
            case 60:
                porder = atoi(optarg);
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
    
    fprintf(stdout, "\nRun SPME\n");
    fprintf(stdout, "----------\n");
    
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

    fprintf(stdout, "num-particles = %d\n", npos);
    fprintf(stdout, "box-size = %f\n", box_size);
    
    int num_rhs = 1;
    double alpha = 1.0;
    double beta = 0.0;
    
    // compute ewald
    MobDebug *ewald = new MobDebug(npos, rdi, box_size,
                                   1.0e-12, MobDebug::EWALD);
    assert(ewald->Init());    
    MobDebug *ewald_real = new MobDebug(npos, rdi, box_size,
                                        1.0e-12, xi, MobDebug::EWALD_REAL);
    assert(ewald_real->Init());
     
    // compute spme
    MobSpme *spme = new MobSpme(npos, rdi, box_size, xi, rmax, dim, porder);
    assert(spme->Init());
    spme->Update(pos, rdi);

    // allocate buffers
    int nm = spme->dim();
    int ldm = nm;
    double *f = (double *)malloc(sizeof(double) * ldm * num_rhs);
    assert(f != NULL);
    RndStream *rnd = new RndStream(1234);
    assert(rnd->Init());
    rnd->Uniform(0.0, 1.0, num_rhs, nm, ldm, f);  
    double *v0 = (double *)malloc(sizeof(double) * ldm * num_rhs);
    double *v0_real = (double *)malloc(sizeof(double) * ldm * num_rhs);
    double *v0_recip = (double *)malloc(sizeof(double) * ldm * num_rhs);
    double *v1 = (double *)malloc(sizeof(double) * ldm * num_rhs);
    double *v1_real = (double *)malloc(sizeof(double) * ldm * num_rhs);
    double *v1_recip = (double *)malloc(sizeof(double) * ldm * num_rhs);
    assert(v0 != NULL && v0_real != NULL && v0_recip != NULL &&
           v1 != NULL && v1_real != NULL && v1_recip != NULL);
    
    // Multiply vectors    
    ewald->MulVector(pos, rdi, num_rhs, alpha, ldm, f, beta, ldm, v0);
    ewald_real->MulVector(pos, rdi, num_rhs, alpha, ldm, f, beta, ldm, v0_real);
    memcpy(v0_recip, v0, sizeof(double) * ldm);
    cblas_daxpy(ldm, -1.0, v0_real, 1, v0_recip, 1);
    
    spme->RealMulVector(num_rhs, alpha, ldm, f, beta, ldm, v1_real);
    spme->RecipMulVector(num_rhs, alpha, ldm, f, beta, ldm, v1_recip);
    memcpy(v1, v1_real, sizeof(double) * ldm);
    cblas_daxpy(ldm, 1.0, v1_recip, 1, v1, 1);
       
    // check results
    MatrixIO *matrix_io = new MatrixIO(); 
    double error;
    matrix_io->Compare(num_rhs, nm, ldm, v0, ldm, v1, &error);
    fprintf(stdout, "\nerror             = %g\n", error);
    matrix_io->Compare(num_rhs, nm, ldm, v0_real, ldm, v1_real, &error);
    fprintf(stdout, "real-space error  = %g\n", error);
    matrix_io->Compare(num_rhs, nm, ldm, v0_recip, ldm, v1_recip, &error);
    fprintf(stdout, "recip-space error = %g\n", error);
    
    // clean up
    delete matrix_io;
    delete parser;
    delete spme;
    delete ewald;
    delete ewald_real;
    delete rnd;
    free(pos);
    free(rdi);
    free(f);
    free(v0);
    free(v0_real);
    free(v0_recip);
    free(v1);
    free(v1_real);
    free(v1_recip);
    free(model_file);
    free(xyz_file);
    
    return 0;
}

