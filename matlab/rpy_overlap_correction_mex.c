#include <math.h>
#include "mex.h"

void rpy_overlap_correction(int np, const double *pos, double L, double *a)
{
    int i, j;
    double posi[4];
    double rvec[4];
    double s2, s;

    int ld = np*3;
    int base;

    // loop over all pairs
    for (i=0; i<np; i++)
    {
        posi[0] = pos[3*i  ];
        posi[1] = pos[3*i+1];
        posi[2] = pos[3*i+2];

        for (j=0; j<np; j++)
        {
            // compute minimum image difference
            rvec[0] = remainder(posi[0] - pos[3*j  ], L);
            rvec[1] = remainder(posi[1] - pos[3*j+1], L);
            rvec[2] = remainder(posi[2] - pos[3*j+2], L);

            s2 = rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2];

            if (s2 < 4 && i != j)
            {
                double t, t1, t2;

                s = sqrt(s2);
                rvec[0] /= s;
                rvec[1] /= s;
                rvec[2] /= s;

                t = 0.09375*s; // 3/32*s
                t1 = (1-3*t) - 0.75/s*(1+2/(3*s2));
                t2 =      t  - 0.75/s*(1-2/s2);

                // entries of a to update are
#define mat(k, l) a[base + ld*k + l]
                base = 3*i*ld + 3*j;
                mat(0,0) = t2*rvec[0]*rvec[0] + t1;
                mat(0,1) = t2*rvec[0]*rvec[1];
                mat(0,2) = t2*rvec[0]*rvec[2];
                mat(1,0) = t2*rvec[1]*rvec[0];
                mat(1,1) = t2*rvec[1]*rvec[1] + t1;
                mat(1,2) = t2*rvec[1]*rvec[2];
                mat(2,0) = t2*rvec[2]*rvec[0];
                mat(2,1) = t2*rvec[2]*rvec[1];
                mat(2,2) = t2*rvec[2]*rvec[2] + t1;
            }
        }
    }
}

// mat = rpy_overlap_correction_mex(pos, L)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int np;
    const double *pos;
    double L;
    double *a;

    np  = mxGetN(prhs[0]);
    pos = mxGetPr(prhs[0]);
    L   = mxGetScalar(prhs[1]);

    // allocate space for output
    plhs[0] = mxCreateDoubleMatrix(np*3, np*3, mxREAL);
    a = mxGetPr(plhs[0]);

    rpy_overlap_correction(np, pos, L, a);
}
