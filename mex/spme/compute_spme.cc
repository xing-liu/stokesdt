#include <stdio.h>
#include <assert.h>
#include <getopt.h>
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
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{    
    if (nrhs < 9) {
        mexErrMsgTxt("compute_spme.m: invalid input.\n");
    }

    // Read prhs[0]
    int npos = mxGetM(prhs[0])/3;
    if (mxGetM(prhs[0])%3 != 0 || npos <= 0) {
        mexErrMsgTxt("compute_spme.m: invalid pos.\n");
    }
    double *pos = mxGetPr(prhs[0]);

    // Read prhs[1]
    if (mxGetM(prhs[1]) != npos) {
        mexErrMsgTxt("compute_spme.m: invalid rdi.\n");
    }
    double *rdi = mxGetPr(prhs[1]);
    
    // Read prhs[2]
    if (mxGetM(prhs[2]) != 1) {
        mexErrMsgTxt("compute_spme.m: invalid L.\n");
    }
    double *data = mxGetPr(prhs[2]);
    double L = data[0];
    if (L <= 0.0) {
        mexErrMsgTxt("compute_spme.m: L is smaller than or equal to 0.0.\n");
    }

    // Read prhs[3]
    if (mxGetM(prhs[3]) != npos * 3) {
        mexErrMsgTxt("compute_spme.m: invalid rdi.\n");
    }
    int num_rhs = mxGetN(prhs[3]);
    double *f = mxGetPr(prhs[3]);

    // Read prhs[4]
    char *mode = mxArrayToString(prhs[4]);
    if (mode == NULL) {
        mexErrMsgTxt("compute_spme.m: invalid input.\n");
    }
    
    // Read prhs[5]
    double xi = -1.0;    
    data = mxGetPr(prhs[5]);    
    xi = data[0]; 
    if (xi <= 0.0) {
        mexErrMsgTxt("compute_spme.m: xi is"
                     " smaller than or equal to 0.0.\n");
    }

    // Read prhs[6]
    double rmax;    
    data = mxGetPr(prhs[6]);    
    rmax = data[0]; 
    if (rmax <= 0.0) {
        mexErrMsgTxt("compute_spme.m: rmax is"
                     " smaller than or equal to 0.0.\n");
    }
    if (2.0 * rmax >= L) {
        mexErrMsgTxt("rmax is larger than or equal to L/2.0");
    }

    // Read prhs[7]
    int dim;    
    data = mxGetPr(prhs[7]);    
    dim = data[0]; 
    if (dim <= 0) {
        mexErrMsgTxt("compute_spme.m: dim is"
                     " smaller than or equal to 0.\n");
    }

    // Read prhs[8]
    int porder;    
    data = mxGetPr(prhs[8]);    
    porder = data[0]; 
    if (porder <= 1) {
        mexErrMsgTxt("compute_spme.m: porder is"
                     " smaller than or equal to 1.\n");
    }
    
    // Create a MobDebug object
    int nm = npos * 3;
    MobSpme *mob = new MobSpme(npos, rdi, L, xi, rmax, dim, porder);
    if (!mob->Init()) {
        mexErrMsgTxt("compute_spme.m: failed to initialize"
                     " the compute engine.\n");    
    }
    mob->Update(pos, rdi);

    // Compute the vector of velocities   
    plhs[0] = mxCreateDoubleMatrix(nm, num_rhs, mxREAL);
    data = mxGetPr(plhs[0]);
    if (strcmp(mode, "full") == 0) {
        mob->MulVector(num_rhs, 1.0, nm, f, 0.0, nm, data);
    } else if (strcmp(mode, "real") == 0) {      
        mob->RealMulVector(num_rhs, 1.0, nm, f, 0.0, nm, data);
    } else if (strcmp(mode, "recip") == 0) {
        mob->RecipMulVector(num_rhs, 1.0, nm, f, 0.0, nm, data);
    } else {
        mexErrMsgTxt("compute_ewald.m: invalid input.\n");
    }
    
    delete mob;
}
