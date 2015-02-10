#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>

#include "force_steric.h"
#include "mex.h"


using namespace stokesdt;

/**
 * prhs[0]: pos(nm x 1)    vector of particle coordinates
 * prhs[1]: rdi(npos x 1)  vector of particle radii
 * prhs[2]: L(1 x 1)       simulation box size
 * prhs[3]: k0(1 x 1)      force constant
 *      
 * plhs[0]: forces(nm x 1) force vector
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{    
    if (nrhs != 4) {
        mexErrMsgTxt("compute_steric_forces.m: invalid input.\n");
    }

    // Read prhs[0]
    int npos = mxGetM(prhs[0])/3;
    if (mxGetM(prhs[0])%3 != 0 || npos <= 0) {
        mexErrMsgTxt("compute_steric_forces.m: invalid pos.\n");
    }
    double *pos = mxGetPr(prhs[0]);

    // Read prhs[1]
    if (mxGetM(prhs[1]) != npos) {
        mexErrMsgTxt("compute_steric_forces.m: invalid rdi.\n");
    }
    double *rdi = mxGetPr(prhs[1]);
    
    // Read prhs[2]
    if (mxGetM(prhs[2]) != 1) {
        mexErrMsgTxt("compute_steric_forces.m: invalid L.\n");
    }
    double *data = mxGetPr(prhs[2]);
    double L = data[0];
    if (L <= 0.0) {
        mexErrMsgTxt("compute_steric_forces.m: L is smaller than or equal to 0.0.\n");
    }

    // Read prhs[3] 
    data = mxGetPr(prhs[3]);    
    double k0 = data[0]; 
    if (k0 <= 0.0) {
        mexErrMsgTxt("compute_steric_forces.m: k0 is"
                     " smaller than or equal to 0.0.\n");
    }

    // Create a StericForce object [what is the overhead of this?]
    // what should be the cutoff
    double r0 = 2.0; // k0 normally 125.
    ForceBase *steric = new StericForce(npos, rdi, L, r0, k0);
    if (!steric->Init()) {
        mexErrMsgTxt("compute_steric_forces.m: Init failed.\n");
    }
    
    // Allocate space and compute
    plhs[0] = mxCreateDoubleMatrix(3*npos, 1, mxREAL);
    steric->Compute(pos, rdi, mxGetPr(plhs[0]));
    
    // Destroy the StericForce object
    delete steric;
}
