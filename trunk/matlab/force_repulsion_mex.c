#include <math.h>
#include "mex.h"

void force_repulsion(int np, const double *pos, double L, double krepulsion, 
    double *forces)
{
    int i, j;
    double posi[4];
    double rvec[4];
    double s2, s, f;

    // initialize forces to zero
    for (i=0; i<3*np; i++)
        forces[i] = 0.;

    // loop over all pairs
    for (i=0; i<np; i++)
    {
        posi[0] = pos[3*i  ];
        posi[1] = pos[3*i+1];
        posi[2] = pos[3*i+2];

        for (j=i+1; j<np; j++)
        {
            // compute minimum image difference
            rvec[0] = remainder(posi[0] - pos[3*j  ], L);
            rvec[1] = remainder(posi[1] - pos[3*j+1], L);
            rvec[2] = remainder(posi[2] - pos[3*j+2], L);

            s2 = rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2];

            if (s2 < 4)
            {
                s = sqrt(s2);
                rvec[0] /= s;
                rvec[1] /= s;
                rvec[2] /= s;
                f = krepulsion*(2.-s);

                forces[3*i  ] +=  f*rvec[0];
                forces[3*i+1] +=  f*rvec[1];
                forces[3*i+2] +=  f*rvec[2];
                forces[3*j  ] += -f*rvec[0];
                forces[3*j+1] += -f*rvec[1];
                forces[3*j+2] += -f*rvec[2];
            }
        }
    }
}

// forces = force_repulsion_mex(pos, L, krepulsion)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int np;
    const double *pos;
    double L;
    double krepulsion;
    double *forces;

    np  = mxGetN(prhs[0]);
    pos = mxGetPr(prhs[0]);
    L   = mxGetScalar(prhs[1]);
    krepulsion = mxGetScalar(prhs[2]);

    // allocate space for output
    plhs[0] = mxCreateDoubleMatrix(3, np, mxREAL);
    forces = mxGetPr(plhs[0]);

    force_repulsion(np, pos, L, krepulsion, forces);
}
