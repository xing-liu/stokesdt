#include <string.h>
#include "mob_spme.h"
#include "mex.h"


using namespace stokesdt;

/**
 * prhs[0]: pos(nm x 1)      the vector of particle coordinates
 * prhs[1]: rdi(npos x 1)    the vector of particle coordinates
 * prhs[2]: L(1 x 1)         the simulation box size
 * plhs[3]: f(nm x num_rhs)  the vector of forces 
 * prhs[4]: mode(1 x 1)      'full', 'real' or 'recip'
 * prhs[5]: xi(1 x 1)        the Ewald parameter
 * prhs[6]: rmax(1 x 1)      the real-space cutoff
 * prhs[7]: dim(1 x 1)       the dimension of FFT grid
 * prhs[8]: porder(1 x 1)    the interpolation order
 *      
 * plhs[1]: v(nm x num_rhs)  the vector of velocities
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    if (nrhs != 9) {
        mexErrMsgTxt("Incorrect number of input arguments.");
    }

    const double *pos    =  mxGetPr(prhs[0]);
    const double *rdi    =  mxGetPr(prhs[1]);
    const double  L      = *mxGetPr(prhs[2]);
    const double *f      =  mxGetPr(prhs[3]);
    const char   *mode   =  mxArrayToString(prhs[4]);
    const double  xi     = *mxGetPr(prhs[5]);
    const double  rmax   = *mxGetPr(prhs[6]);
    const int     dim    = (int ) *mxGetPr(prhs[7]);
    const int     porder = (int ) *mxGetPr(prhs[8]);

    const int npos = mxGetN(prhs[0]);
    const int num_rhs = mxGetN(prhs[3]);

    if (mxGetM(prhs[0]) != 3) {
        mexErrMsgTxt("Invalid argument (pos must have leading dimension 3).");
    }

    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != npos) {
        mexErrMsgTxt("Different number of particles for rdi and pos.");
    }

    if (L <= 0. || xi <= 0. || rmax <= 0 || dim < 1) {
        mexErrMsgTxt("Scalar arguments must be positive.");
    }

    if (porder <= 1) {
        mexErrMsgTxt("porder must be greater than 1.");
    }

    if (2.*rmax >= L) {
        mexErrMsgTxt("rmax must be less than L/2.");
    }

    if (mode == NULL) {
        mexErrMsgTxt("invalid argument: mode.");
    }

    // Create a MobDebug object
    MobSpme *mob = new MobSpme(npos, rdi, L, xi, rmax, dim, porder);
    if (!mob->Init()) {
        mexErrMsgTxt("Init failed.");
    }
    mob->Update(pos, rdi);

    // Compute the vector of velocities   
    int nm = 3*npos;
    plhs[0] = mxCreateDoubleMatrix(nm, num_rhs, mxREAL);
    double *data = mxGetPr(plhs[0]);
    if (strcmp(mode, "full") == 0) {
        mob->MulVector(num_rhs, 1.0, nm, f, 0.0, nm, data);
    } else if (strcmp(mode, "real") == 0) {      
        mob->RealMulVector(num_rhs, 1.0, nm, f, 0.0, nm, data);
    } else if (strcmp(mode, "recip") == 0) {
        mob->RecipMulVector(num_rhs, 1.0, nm, f, 0.0, nm, data);
    } else {
        mexErrMsgTxt("mode must be full, real, or recip.");
    }
    
    delete mob;
}
