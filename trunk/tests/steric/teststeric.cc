#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <cstring>

#include "stokes_force.h"
#include "stokes_util.h"


using namespace stokesdt;

const struct option long_options[] = {
    {"model",   1, NULL, 10},
    {"xyz",     1, NULL, 20},
    {"cutoff",  1, NULL, 30},
    {"k0",      1, NULL, 40},
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
    fprintf(stderr, "\t--cutoff             Steric cutoff\n");    
    fprintf(stderr, "\t--k0                 Steric constant\n");
}


static int DirectSteric(const double box_size,
                        const int npos,
                        const double *pos,
                        const double *rdi,
                        const double cutoff,
                        const double k0,
                        double *f)
{
    int nnz = 0;
    memset(f, 0, sizeof(double) * npos * 3);

    #pragma omp parallel for reduction(+:nnz)
    for (int i = 0; i < npos; i++) {
        double xi = pos[i * 3];
        double yi = pos[i * 3 + 1];
        double zi = pos[i * 3 + 2];
        double ai = rdi[i];
        for (int j = 0; j < npos; j++) {
            if (i == j)
                continue;
            double xj = pos[j * 3];
            double yj = pos[j * 3 + 1];
            double zj = pos[j * 3 + 2];
            double aj = rdi[j];
            double rx = drem(xi - xj, box_size);
            double ry = drem(yi - yj, box_size);
            double rz = drem(zi - zj, box_size);
            double rr = rx * rx + ry * ry + rz * rz;
            double r = sqrt(rr);
            r = (r == 0.0 ? 1.0e-10 : r);
            double ex = rx / r;
            double ey = ry / r;
            double ez = rz / r;
            double s = 2.0 * r / (ai + aj);
            if (s - cutoff >= 0.0)
                continue;
            double force = k0 * ((ai + aj) * cutoff / 2.0 - r);
            f[3 * i + 0] +=  force * ex;
            f[3 * i + 1] +=  force * ey;
            f[3 * i + 2] +=  force * ez;
            nnz++;
        }
    }
     
    return nnz + npos;
}


int main(int argc, char **argv)
{
    char *model_file = NULL;
    char *xyz_file = NULL;
    double cutoff = 2.0;
    double k0 = 125.0;   
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
                cutoff = atof(optarg);
                break;                
            case 40:
                k0 = atof(optarg);
                break;                
            case ':':
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                return -1;
            case '?':
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
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
    
    fprintf(stdout, "\nTest steric force\n");
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

    // compute steric
    int nf = npos * 3;
    double *f = (double *)malloc(sizeof(double) * nf);
    assert(f != NULL);
    StericForce *steric = new StericForce(npos, rdi, box_size, cutoff, k0);
    assert(steric->Init());
    steric->Compute(pos, rdi, f);

    // check result
    double *f0 = (double *)malloc(sizeof(double) * nf);
    assert(f0 != NULL);
    DirectSteric(box_size, npos, pos, rdi, cutoff, k0, f0);
    MatrixIO *matrix_io = new MatrixIO();
    double error;
    matrix_io->Compare(1, nf, nf, f0, nf, f, &error);
    fprintf(stdout, "error = %le\n", error);
    if (error < 1.0e-6) {
        fprintf(stdout, "PASS.\n");
    } else {
        fprintf(stdout, "FAILED.\n\n");
    }
    
    free(f);
    free(f0);
    delete steric;
    delete parser;
    delete matrix_io;
    free(model_file);
    free(xyz_file);
        
    return 0;
}
