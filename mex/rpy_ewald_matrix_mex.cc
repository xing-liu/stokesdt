#include <string.h>
#include "mob_debug.h"
#include "mex.h"


using namespace stokesdt;

// mat = rpy_ewald_matrix_mex(pos, rdi, L, tol, mode, xi)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    if (nrhs != 5 && nrhs != 6) {
        mexErrMsgTxt("Needs 5 or 6 input arguments.");
    }

    const double *pos    =  mxGetPr(prhs[0]);
    const double *rdi    =  mxGetPr(prhs[1]);
    const double  L      = *mxGetPr(prhs[2]);
    const double  tol    = *mxGetPr(prhs[3]);
    const char   *mode   =  mxArrayToString(prhs[4]);
          double  xi; // will be set below

    const int npos = mxGetN(prhs[0]);

    if (mxGetM(prhs[0]) != 3) {
        mexErrMsgTxt("Invalid argument (pos must have leading dimension 3).");
    }

    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != npos) {
        mexErrMsgTxt("Different number of particles for rdi and pos.");
    }

    if (L <= 0. || tol <= 0.) {
        mexErrMsgTxt("Scalar arguments must be positive.");
    }

    if (mode == NULL) {
        mexErrMsgTxt("Invalid argument: mode.");
    }

    MobDebug::MobDebugType mob_type;
    if (strcmp(mode, "full") == 0) {
        mob_type = MobDebug::EWALD;
    } else if (strcmp(mode, "real") == 0) {      
        mob_type = MobDebug::EWALD_REAL;
    } else if (strcmp(mode, "recip") == 0) {
        mob_type = MobDebug::EWALD_RECIP;
    } else {
        mexErrMsgTxt("mode must be full, real, or recip.");
    }

    // set xi
    xi = (nrhs == 6) ? (*mxGetPr(prhs[5])) : (pow(10.0, 1.0/6.0)*sqrt(M_PI)/L);

    MobDebug *mob = new MobDebug(npos, rdi, L, tol, xi, mob_type);
    if (!mob->Init()) {
        mexErrMsgTxt("Init failed.");
    }
    
    plhs[0] = mxCreateDoubleMatrix(3*npos, 3*npos, mxREAL);
    mob->GetMatrix(pos, rdi, 3*npos, mxGetPr(plhs[0]));
    
    delete mob;
}
