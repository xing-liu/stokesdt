#include <stdio.h>
#include <math.h>
#include <omp.h>

#define WALLCLOCK(time) do {                                 \
      unsigned long val;                                       \
      volatile unsigned int a, d;                              \
      __asm__ __volatile__("rdtsc" : "=a" (a), "=d" (d) : );   \
      val = ((unsigned long) a)|(((unsigned long)d)<<32);      \
      (time) = val / 2400000000.;                              \
    } while(0)

// this should be vectorized
void scalar_rpy_ewald_real(double r, double xi, double *m11, double *m12)
{
    double a = 1.;
    double a3 = 1.;

    double xi2 = xi*xi;
    double xi3 = xi2*xi;
    double xi5 = xi3*xi2;
    double xi7 = xi5*xi2;

    double r2 = r*r;
    double r4 = r2*r2;
    double ri = 1./r;
    double ri2 = ri*ri;
    double ri3 = ri*ri2;

    double erfc_xi_r = erfc(xi*r);
    double pi_exp = 1./sqrt(M_PI) * exp(-xi2*r2);

    *m11 = (0.75*a*ri + 0.5*a3*ri3)*erfc_xi_r + ( 4*xi7*a3*r4 + 3*xi3*a*r2 - 20*xi5*a3*r2 - 4.5*xi*a + 14*xi3*a3 +   xi*a3*ri2)*pi_exp;
    *m12 = (0.75*a*ri - 1.5*a3*ri3)*erfc_xi_r + (-4*xi7*a3*r4 - 3*xi3*a*r2 + 16*xi5*a3*r2 + 1.5*xi*a -  2*xi3*a3 - 3*xi*a3*ri2)*pi_exp;
}

void scalar_rpy_ewald_recip(double k, double xi, double *m2)
{
    double a = 1.;
    double a3 = 1.;

    double k2 = k*k;
    double xii2k2 = k2/(xi*xi);

    *m2 = (a-a3*k2/3.) * (1. + 0.25*xii2k2 + 0.125*xii2k2*xii2k2) * 6.*M_PI/k2 * exp(-0.25*xii2k2);
}

// note: positions must be wrapped inside the box [0,L]
int rpy_ewald(int np, double *a, 
    const double *pos, double L, double xi, int nr, int nk)
{
    double posi[4];
    double posj[4];
    double rvec[4];
    double rvec0[4];
    double temp[6];

    double m11, m12, m2;
    double eye3_coef;
    double r2, r;

    int x, y, z;
    int i, j;

    double *ap0, *ap;

    int vsize = ((2*nk+1)*(2*nk+1)*(2*nk+1) - 1) / 2;
#define VSIZE ((2*6+1)*(2*6+1)*(2*6+1) - 1) / 2
    double k_array[VSIZE];
    double m2_array[VSIZE];
    double kvec_array[3*VSIZE];
    int ind;

    double kvec[4];
    double k;
    double t;

    double vinv = 1./(L*L*L);

    double time0, time1;
    double time0_real, time1_real;
    double time0_recip, time1_recip;
    WALLCLOCK(time0);

    ind = 0;
    for (x=-nk; x<=nk; x++)
    for (y=-nk; y<=nk; y++)
    for (z=-nk; z<=nk; z++)
    {
        k_array[ind] = 2.*M_PI/L*sqrt((double)(x*x + y*y + z*z));
        scalar_rpy_ewald_recip(k_array[ind], xi, &m2_array[ind]);
        kvec_array[3*ind  ] = 2.*M_PI/L*x;
        kvec_array[3*ind+1] = 2.*M_PI/L*y;
        kvec_array[3*ind+2] = 2.*M_PI/L*z;
        ind++;
        if (ind >= vsize)
            goto out;
    }
out:
    // *************************************************************************
    // real-space sum

    WALLCLOCK(time0_real);
#pragma omp parallel num_threads(20) private(i,j,posi,posj,temp,rvec,rvec0,x,y,z,r2,r,eye3_coef,m11,m12,ap0,ap,   kvec,k,t,m2,ind)
  {
#pragma omp for schedule(dynamic)
    for (i=0; i<np; i++)
    {
        posi[0] = pos[3*i];
        posi[1] = pos[3*i+1];
        posi[2] = pos[3*i+2];

        for (j=i; j<np; j++)
        {
            temp[0] = 0.;
            temp[1] = 0.;  temp[3] = 0.;
            temp[2] = 0.;  temp[4] = 0.;  temp[5] = 0.;
            eye3_coef = 0.;

            posj[0] = pos[3*j];
            posj[1] = pos[3*j+1];
            posj[2] = pos[3*j+2];

            rvec0[0] = posi[0] - posj[0];
            rvec0[1] = posi[1] - posj[1];
            rvec0[2] = posi[2] - posj[2];

            for (x=-nr; x<=nr; x++)
            for (y=-nr; y<=nr; y++)
            for (z=-nr; z<=nr; z++)
            {
                rvec[0] = rvec0[0] + x*L;
                rvec[1] = rvec0[1] + y*L;
                rvec[2] = rvec0[2] + z*L;

                // compute norm
                r2 = rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2];
                r  = sqrt(r2);

                if (r2 == 0.)
                    continue;

                // rvec = rvec/r
                rvec[0] /= r;
                rvec[1] /= r;
                rvec[2] /= r;

                // compute scalars for this value of r
                scalar_rpy_ewald_real(r, xi, &m11, &m12);

                // we can sum the coeff for m11 outside 
                eye3_coef += m11;

                temp[0] += m12 * rvec[0] * rvec[0];
                temp[1] += m12 * rvec[0] * rvec[1];
                temp[2] += m12 * rvec[0] * rvec[2];
                temp[3] += m12 * rvec[1] * rvec[1];
                temp[4] += m12 * rvec[1] * rvec[2];
                temp[5] += m12 * rvec[2] * rvec[2];
            }

            // add contribution to eye3 term
            temp[0] += eye3_coef;
            temp[3] += eye3_coef;
            temp[5] += eye3_coef;

            if (i == j)
            {
               temp[0] *= 0.5;
               temp[1] *= 0.5;
               temp[2] *= 0.5;
               temp[3] *= 0.5;
               temp[4] *= 0.5;
               temp[5] *= 0.5;
            }

            // sum into global matrix (only lower-triangular part)
            // Use matlab to add transpose

             ap0  = &a[np*3*3*i + 3*j];
             ap   = ap0;
            *ap++ = temp[0];
            *ap++ = temp[1];
            *ap   = temp[2];
             ap   = ap0+np*3;
            *ap++ = temp[1];
            *ap++ = temp[3];
            *ap   = temp[4];
             ap   = ap0+np*3+np*3;
            *ap++ = temp[2];
            *ap++ = temp[4];
            *ap   = temp[5];
        }
    }
    WALLCLOCK(time1_real);

    // *************************************************************************
    // compute and save coefficients for reciprocal-space sum
    // Due to symmetry, only need half of the grid points

#if 0
    ind = 0;
    for (x=-nk; x<=nk; x++)
    for (y=-nk; y<=nk; y++)
    for (z=-nk; z<=nk; z++)
    {
        k_array[ind] = 2.*M_PI/L*sqrt((double)(x*x + y*y + z*z));
        scalar_rpy_ewald_recip(k_array[ind], xi, &m2_array[ind]);
        kvec_array[3*ind  ] = 2.*M_PI/L*x;
        kvec_array[3*ind+1] = 2.*M_PI/L*y;
        kvec_array[3*ind+2] = 2.*M_PI/L*z;
        ind++;
        if (ind >= vsize)
            goto out;
    }
out:
#endif

    // *************************************************************************
    // reciprocal-space sum

    WALLCLOCK(time0_recip);
#pragma omp for schedule(dynamic)
    for (i=0; i<np; i++)
    {
        posi[0] = pos[3*i];
        posi[1] = pos[3*i+1];
        posi[2] = pos[3*i+2];

        for (j=i; j<np; j++)
        {
            rvec[0] = posi[0] - pos[3*j];
            rvec[1] = posi[1] - pos[3*j+1];
            rvec[2] = posi[2] - pos[3*j+2];

            temp[0] = 0.;
            temp[1] = 0.;  temp[3] = 0.;
            temp[2] = 0.;  temp[4] = 0.;  temp[5] = 0.;

            for (ind=0; ind<vsize; ind++)
            {
                k = k_array[ind];
                m2 = m2_array[ind];
                kvec[0] = kvec_array[3*ind  ];
                kvec[1] = kvec_array[3*ind+1];
                kvec[2] = kvec_array[3*ind+2];

                t = 2.*vinv*m2*cos(
                    kvec[0]*rvec[0] + kvec[1]*rvec[1] + kvec[2]*rvec[2]);

                kvec[0] /= k;
                kvec[1] /= k;
                kvec[2] /= k;

                temp[0] += t * (1. - kvec[0]*kvec[0]);
                temp[1] += t *     - kvec[0]*kvec[1];
                temp[2] += t *     - kvec[0]*kvec[2];
                temp[3] += t * (1. - kvec[1]*kvec[1]);
                temp[4] += t *     - kvec[1]*kvec[2];
                temp[5] += t * (1. - kvec[2]*kvec[2]);
            }

            // sum into matrix
            if (i == j)
            {
                temp[0] *= 0.5;
                temp[1] *= 0.5;
                temp[2] *= 0.5;
                temp[3] *= 0.5;
                temp[4] *= 0.5;
                temp[5] *= 0.5;
            }

            // sum with existing values
             ap0   = &a[np*3*3*i + 3*j];
             ap    = ap0;
            *ap++ += temp[0];
            *ap++ += temp[1];
            *ap   += temp[2];
             ap    = ap0+np*3;
            *ap++ += temp[1];
            *ap++ += temp[3];
            *ap   += temp[4];
             ap    = ap0+np*3+np*3;
            *ap++ += temp[2];
            *ap++ += temp[4];
            *ap   += temp[5];
        }
    }
    WALLCLOCK(time1_recip);
    
  }

    // *************************************************************************
    // self-part

    double radius = 1;
    double radius3 = radius*radius*radius;

    t = 1. - (6.*radius - 40./3.*xi*xi*radius3)*xi/sqrt(M_PI);
    t *= 0.5;
    for (i=0; i<np*3; i++)
        a[i*np*3+i] += t;

    WALLCLOCK(time1);
    //printf("%d %d wallclock real recip total: %f %f %f\n", 
#if 0
    printf("%d %d & %f %f %f & ", 
        nr, nk,
        time1_real-time0_real,
        time1_recip-time0_recip,
        time1-time0);
#endif

    return 0;
}

#if 1
#include "mex.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    int np;
    double *a;
    const double *pos;
    double L;
    double xi;
    int nr;
    int nk;

    // arguments should have been checked in m-file
    np  = mxGetN(prhs[0]);
    pos = mxGetPr(prhs[0]);
    L   = mxGetScalar(prhs[1]);
    xi  = mxGetScalar(prhs[2]);
    nr  = (int) mxGetScalar(prhs[3]);
    nk  = (int) mxGetScalar(prhs[4]);

    // allocate space for output
    plhs[0] = mxCreateDoubleMatrix(np*3, np*3, mxREAL);
    a = mxGetPr(plhs[0]);

    rpy_ewald(np, a, pos, L, xi, nr, nk);
}

#else
int main()
{
    int np = 3;
    double pos[] = {0., 0., 0.,  0., 0., 4.,  0., 4., 0.};
    double phi = 0.5;
    double L = pow(4./3.*M_PI*np/phi, 1./3.);
    double xi = sqrt(M_PI)/L;
    int nr = 4;
    int nk = 4;

    double a[3*3*3*3];

    rpy_ewald(np, a, pos, L, xi, nr, nk);

    return 0;
}
#endif

