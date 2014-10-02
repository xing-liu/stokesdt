#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>

#include "mob_debug.h"
#include "mex.h"


using namespace stokesdt;

/**
 * prhs[0]: pos(nm x 1)    the vector of particle coordinates
 * prhs[1]: rdi(npos x 1)  the vector of particle coordinates
 * prhs[2]: L(1 x 1)       the simulation box size
 * prhs[3]: tol(1 x 1)     the requested tolerance of Ewald errors
 * prhs[4]: mode(1 x 1)    'full', 'real' or 'recip'
 * prhs[5]: xi(1 x 1)      the Ewald parameter (optional);
 *                         if not specified, the optimal xi will be
 *                         automatically chosen
 *      
 * plhs[1]: mat(nm by nm)  the array of the mobility matrix
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{    
    if (nrhs < 5) {
        mexErrMsgTxt("compute_ewald.m: invalid input.\n");
    }

    // Read prhs[0]
    int npos = mxGetM(prhs[0])/3;
    if (mxGetM(prhs[0])%3 != 0 || npos <= 0) {
        mexErrMsgTxt("compute_ewald.m: invalid pos.\n");
    }
    double *pos = mxGetPr(prhs[0]);

    // Read prhs[1]
    if (mxGetM(prhs[1]) != npos) {
        mexErrMsgTxt("compute_ewald.m: invalid rdi.\n");
    }
    double *rdi = mxGetPr(prhs[1]);
    
    // Read prhs[2]
    if (mxGetM(prhs[2]) != 1) {
        mexErrMsgTxt("compute_ewald.m: invalid L.\n");
    }
    double *data = mxGetPr(prhs[2]);
    double L = data[0];
    if (L <= 0.0) {
        mexErrMsgTxt("compute_ewald.m: L is smaller than or equal to 0.0.\n");
    }

    // Read prhs[3]
    if (mxGetM(prhs[3]) != 1) {
        mexErrMsgTxt("compute_ewald.m: invalid input.\n");
    }
    data = mxGetPr(prhs[3]);    
    double tol = data[0]; 
    if (tol <= 0.0) {
        mexErrMsgTxt("compute_ewald.m: tol is smaller than or equal to 0.0.\n");
    }

    // Read prhs[4]
    char *mode = mxArrayToString(prhs[4]);
    if (mode == NULL) {
        mexErrMsgTxt("compute_ewald.m: invalid input.\n");
    }
    MobDebug::MobDebugType mob_type;
    if (strcmp(mode, "full") == 0) {
        mob_type = MobDebug::EWALD;
    } else if (strcmp(mode, "real") == 0) {      
        mob_type = MobDebug::EWALD_REAL;
    } else if (strcmp(mode, "recip") == 0) {
        mob_type = MobDebug::EWALD_RECIP;
    } else {
        mexErrMsgTxt("compute_ewald.m: invalid input.\n");
    }

    // Read prhs[5]
    double xi = -1.0;
    if (nrhs == 6) {
        data = mxGetPr(prhs[5]);    
        xi = data[0]; 
        if (xi <= 0.0) {
            mexErrMsgTxt("compute_ewald.m: xi is"
                         " smaller than or equal to 0.0.\n");
        }    
    }
    
    // Create a MobDebug object
    int nm = npos * 3;
    MobDebug *mob = new MobDebug(npos, rdi, L, tol, mob_type);
    if (!mob->Init()) {
        mexErrMsgTxt("compute_ewald.m: failed to initialize"
                     " the compute engine.\n");    
    }
    
    // Write plhs[0]
    plhs[0] = mxCreateDoubleMatrix(nm, nm, mxREAL);
    data = mxGetPr(plhs[0]);
    mob->GetMatrix(pos, rdi, nm, data);
    
    delete mob;
}
