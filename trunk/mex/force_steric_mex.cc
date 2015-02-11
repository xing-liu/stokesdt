#include "force_steric.h"
#include "mex.h"


using namespace stokesdt;

// f = forces_steric_mex(pos, rdi, L, r0, k0)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    if (nrhs != 5) {
        mexErrMsgTxt("Incorrect number of input arguments.");
    }

    double *pos =  mxGetPr(prhs[0]);
    double *rdi =  mxGetPr(prhs[1]);
    double  L   = *mxGetPr(prhs[2]);
    double  r0  = *mxGetPr(prhs[3]);
    double  k0  = *mxGetPr(prhs[4]);

    int npos = mxGetN(prhs[0]);

    if (mxGetM(prhs[0]) != 3) {
        mexErrMsgTxt("Invalid argument (pos must have leading dimension 3).");
    }

    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != npos) {
        mexErrMsgTxt("Different number of particles for rdi and pos.");
    }

    if (L <= 0. || r0 <= 0. || k0 <= 0.) {
        mexErrMsgTxt("Scalar arguments must be positive.");
    }

    ForceBase *steric = new StericForce(npos, rdi, L, r0, k0);
    if (!steric->Init()) {
        mexErrMsgTxt("Init failed.");
    }
    
    plhs[0] = mxCreateDoubleMatrix(3, npos, mxREAL);
    steric->Compute(pos, rdi, mxGetPr(plhs[0]));
    
    delete steric;
}
