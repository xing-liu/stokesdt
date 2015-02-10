#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>

#include "pair_r.h"
#include "pair_s.h"
#include "mex.h"


using namespace stokesdt;

/**
 * prhs[0]: pos(nm x 1)    the vector of particle coordinates
 * prhs[1]: rdi(npos x 1)  the vector of particle coordinates
 * prhs[2]: L(1 x 1)       the simulation box size
 * prhs[3]: cutoff(1 x 1)  the cutoff
 * prhs[4]: mode(1 x 1)    'R' (absulote distance) or 'S' (normalized distance)
 *      
 * plhs[0]: indi             the vector of source particle identities
 * plhs[1]: indj             the vector of destination particle identities
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{    
    if (nrhs != 5) {
        mexErrMsgTxt("compute_pairlist.m: invalid input.\n");
    }

    // Read prhs[0]
    int npos = mxGetM(prhs[0])/3;
    if (mxGetM(prhs[0])%3 != 0 || npos <= 0) {
        mexErrMsgTxt("compute_pairlist.m: invalid pos.\n");
    }
    double *pos = mxGetPr(prhs[0]);

    // Read prhs[1]
    if (mxGetM(prhs[1]) != npos) {
        mexErrMsgTxt("compute_pairlist.m: invalid rdi.\n");
    }
    double *rdi = mxGetPr(prhs[1]);
    
    // Read prhs[2]
    if (mxGetM(prhs[2]) != 1) {
        mexErrMsgTxt("compute_pairlist.m: invalid L.\n");
    }
    double *data = mxGetPr(prhs[2]);
    double L = data[0];
    if (L <= 0.0) {
        mexErrMsgTxt("compute_pairlist.m: L is smaller than or equal to 0.0.\n");
    }

    // Read prhs[3] 
    data = mxGetPr(prhs[3]);    
    double cutoff = data[0]; 
    if (cutoff <= 0.0) {
        mexErrMsgTxt("compute_spme.m: cutoff is"
                     " smaller than or equal to 0.0.\n");
    }
    if (2.0 * cutoff >= L) {
        mexErrMsgTxt("cutoff is larger than or equal to L/2.0");
    }

    // Read prhs[4]
    char *mode = mxArrayToString(prhs[4]);
    if (mode == NULL) {
        mexErrMsgTxt("compute_pairlist.m: invalid input.\n");
    }
    
    // Create Pair list
    PairListBase *pair;
    if (strcmp(mode, "R") == 0) {
        pair = new PairListR(npos, rdi, L, cutoff);
    } else if (strcmp(mode, "S") == 0) {      
        pair = new  PairListS(npos, rdi, L, cutoff);
    } else {        
        mexErrMsgTxt("compute_pairlist.m: invalid input.\n");
    }

    if (pair->Init() <= 0) {
        mexErrMsgTxt("compute_pairlist.m: failed to initialize"
                     " the compute engine.\n");    
    }
    int nnz = pair->Build(pos);
      
    // Write plhs[0], plhs[1]
#if 0
    plhs[0] = mxCreateNumericMatrix(nnz, 1, mxINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(nnz, 1, mxINT32_CLASS, mxREAL);
    int *id1 = (int *)mxGetPr(plhs[0]);
    int *id2 = (int *)mxGetPr(plhs[1]);
    std::vector<int> rowptr (npos + 1);
    pair->GetPairs(&rowptr[0], id2);
    for (int i = 0; i < npos; i++) {
        int start = rowptr[i];
        int end = rowptr[i + 1];
        for (int j = start; j < end; j++) {
            id1[j] = i;
        }
    }
#else
    plhs[0] = mxCreateNumericMatrix(npos+1, 1, mxINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(nnz, 1, mxINT32_CLASS, mxREAL);
    int *id1 = (int *)mxGetPr(plhs[0]);
    int *id2 = (int *)mxGetPr(plhs[1]);
    pair->GetPairs(id1, id2);
#endif
    
    delete pair;
}
