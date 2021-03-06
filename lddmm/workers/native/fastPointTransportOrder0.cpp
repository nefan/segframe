/*
 *
 * This file is part of segframe.
 *
 * Copyright (C) 2012, Stefan Sommer (sommer@diku.dk)
 * https://github.com/nefan/segframe.git
 *
 * segframe is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * segframe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with segframe.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "common.h"
#include "sla.h"

using namespace sla;

/* The computational routine */
void fastPointTransportOrder0(scalar *dgrid, scalar t, scalar *gridt, scalar *rhot, int L, int R, int CSP, scalar *scales2, scalar *scaleweight2, int Ngrid)
{
    /* dgrid already initialized to zeros */
    int i, l, s;
    const int dim = 3;

#pragma omp parallel for schedule(static) shared(dgrid,t,gridt,rhot,L,R,dim,CSP,scales2,scaleweight2) private(i,l,s)
    for (i=0; i<Ngrid; i++) { /* grid point */
        Vector3<scalar> xi(gridt[dim*i+0],gridt[dim*i+1],gridt[dim*i+2]);

        for (l=0; l<L; l++) { /* particle */
            Vector3<scalar> xl(rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)],rhot[INDRHOX(l,2)]);

            Vector3<scalar> ximxl = xi-xl;
            scalar v = dot<scalar>(ximxl,ximxl);

            for (s=0; s<R; s++) { /* scale */
                Vector3<scalar> al(rhot[INDRHOP(l,0,s)],rhot[INDRHOP(l,1,s)],rhot[INDRHOP(l,2,s)]);
                scalar r2 = scales2[s];
                scalar sw2 = scaleweight2[s];
                scalar kbase = exp(-v/r2);
                scalar ks = KKs(kbase,sw2);

                Vector3<scalar> vt1 = ks*al;
                dgrid[dim*i+0] += vt1[0];
                dgrid[dim*i+1] += vt1[1];
                dgrid[dim*i+2] += vt1[2];
            }
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double t;
    double *gridt;
    double *rhot;
    int L;
    int R;
    int dim;
    double *scales2;
    double *scaleweight2;
    double *dgrid;

    /* check for proper number of arguments */
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("gradC:nrhs","8 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("gradC:nlhs","One output required.");
    }
    
    /* get the values of the non-array inputs  */
    t = mxGetScalar(prhs[0]);
    L = (double)(mxGetScalar(prhs[3]));
    R = (double)(mxGetScalar(prhs[4]));
    dim = (double)(mxGetScalar(prhs[5]));
    mexAssert(dim == 3);

    int CSP = dim*(1+R);

    /* get the array inputs */
    gridt = mxGetPr(prhs[1]);
    int n= mxGetM(prhs[1]);
    int Ngrid = n/3;
    rhot = mxGetPr(prhs[2]);
    if(CSP*L != mxGetM(prhs[2])) {
        mexErrMsgIdAndTxt("gradC:nlhs","rhot dimension mismatch (got %d, expected %d).",mxGetM(prhs[2]),CSP*L);
    }
    scales2 = mxGetPr(prhs[6]);
    if(R != mxGetN(prhs[6])) {
        mexErrMsgIdAndTxt("gradC:nlhs","scales2 dimension mismatch (got %d, expected %d).",mxGetN(prhs[6]),R);
    }
    scaleweight2 = mxGetPr(prhs[7]);
    if(R != mxGetN(prhs[7])) {
        mexErrMsgIdAndTxt("gradC:nlhs","scaleweight2 dimension mismatch (got %d, expected %d).",mxGetN(prhs[7]),R);
    }

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    dgrid = mxGetPr(plhs[0]);

    /* debug */
    /*mexPrintf("%f %f %f %d %d %d %d %f %f %d %d %d\n",t,yt[0],rhot[0],L,R,CSP,CSD,scales2[0],scaleweight2[0],GL,GS,GI);*/

    /* call the computational routine */
    fastPointTransportOrder0(dgrid,t,gridt,rhot,L,R,CSP,scales2,scaleweight2,Ngrid);
}
