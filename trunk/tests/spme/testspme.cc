#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <cstring>

#include "stokes_mob.h"
#include "stokes_util.h"


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
    
    fprintf(stdout, "\nTest SPME\n");
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

    int num_rhs = 3;
    double alpha = 2.1;
    double beta = 1.1;
    
    // compute ewald
    MobEwald *mob_ewald = new MobEwald(npos, rdi, box_size, 1.0e-12);
    assert(mob_ewald->Init());
    mob_ewald->Update(pos, rdi);
     
    // compute spme
    MobSpme *mob_spme = new MobSpme(npos, rdi, box_size, xi, rmax, dim, porder);
    assert(mob_spme->Init());
    mob_spme->Update(pos, rdi);

    // allocate buffers
    int nm = mob_spme->dim();
    int ldm = nm;
    double *f = (double *)malloc(sizeof(double) * ldm * num_rhs);
    assert(f != NULL);
    RndStream *rnd = new RndStream(1234);
    assert(rnd->Init());
    rnd->Uniform(0.0, 1.0, num_rhs, nm, ldm, f);  
    double *v = (double *)malloc(sizeof(double) * ldm * num_rhs);
    double *v_spme = (double *)malloc(sizeof(double) * ldm * num_rhs);
    assert(v != NULL && v_spme != NULL);
    memcpy(v, f, sizeof(double) * ldm * num_rhs);
    memcpy(v_spme, f, sizeof(double) * ldm * num_rhs);
    
    // Multiply vectors   
    mob_ewald->MulVector(num_rhs, alpha, ldm, f, beta, ldm, v);
    mob_spme->MulVector(num_rhs, alpha, ldm, f, beta, ldm, v_spme);
   
    // check results
    MatrixIO *matrix_io = new MatrixIO(); 
    double error;
    matrix_io->Compare(num_rhs, nm, ldm, v, ldm, v_spme, &error);
    fprintf(stdout, "error = %le\n", error);

    if (error < 1.0e-6) {
        fprintf(stdout, "PASS.\n");
    } else {
        fprintf(stdout, "FAILED.\n\n");
    }
    
    // clean up
    delete matrix_io;
    delete parser;
    delete mob_spme;
    delete mob_ewald;
    delete rnd;
    free(pos);
    free(rdi);
    free(f);
    free(v);
    free(v_spme);
    free(model_file);
    free(xyz_file);
    
    return 0;
}
