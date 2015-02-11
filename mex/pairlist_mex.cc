#include <string.h>
#include "pair_r.h"
#include "pair_s.h"
#include "mex.h"


using namespace stokesdt;

// [ia, ja] = pairlist_mex(pos, rdi, L, cutoff, mode)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    if (nrhs != 5) {
        mexErrMsgTxt("Incorrect number of input arguments.");
    }

    const double *pos    =  mxGetPr(prhs[0]);
    const double *rdi    =  mxGetPr(prhs[1]);
    const double  L      = *mxGetPr(prhs[2]);
    const double  cutoff = *mxGetPr(prhs[3]);
    const char   *mode   =  mxArrayToString(prhs[4]);

    const int npos = mxGetN(prhs[0]);

    if (mxGetM(prhs[0]) != 3) {
        mexErrMsgTxt("Invalid argument (pos must have leading dimension 3).");
    }

    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != npos) {
        mexErrMsgTxt("Different number of particles for rdi and pos.");
    }

    if (L <= 0. || cutoff <= 0.) {
        mexErrMsgTxt("Scalar arguments must be positive.");
    }

    if (2.*cutoff >= L) {
        mexErrMsgTxt("cutoff must be less than L/2.");
    }

    if (mode == NULL) {
        mexErrMsgTxt("invalid argument: mode.");
    }
    
    // Create Pair list
    PairListBase *pair;
    if (strcmp(mode, "R") == 0) {
        pair = new PairListR(npos, rdi, L, cutoff);
    } else if (strcmp(mode, "S") == 0) {      
        pair = new PairListS(npos, rdi, L, cutoff);
    } else {        
        mexErrMsgTxt("mode must be R or S.");
    }

    if (pair->Init() <= 0) {
        mexErrMsgTxt("Init failed.");
    }
    int nnz = pair->Build(pos);
      
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
    int *ia = (int *)mxGetPr(plhs[0]);
    int *ja = (int *)mxGetPr(plhs[1]);
    pair->GetPairs(ia, ja);
#endif
    
    delete pair;
}
