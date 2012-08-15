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
void fastPointPathOrder0(scalar *drho, scalar t, scalar *rhot, int L, int R, int CSP, scalar *scales2, scalar *scaleweight2)
{
    /* drho already initialized to zeros */
    int i, l, si, sl;
    const int dim = 3;

    for (i=0; i<L; i++) { /* particle */
        Vector3<scalar> xi(rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)],rhot[INDRHOX(i,2)]);

        for (l=0; l<L; l++) { /* particle */
            Vector3<scalar> xl(rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)],rhot[INDRHOX(l,2)]);

            Vector3<scalar> ximxl = xi-xl;
            scalar v = dot<scalar>(ximxl,ximxl);

            for (sl=0; sl<R; sl++) { /* scale */
                Vector3<scalar> alsl(rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)],rhot[INDRHOP(l,2,sl)]);
                scalar rsl2 = scales2[sl];
                scalar swsl2 = scaleweight2[sl];
                scalar kbasesl = exp(-v/rsl2);
                scalar ksl = KKs(kbasesl,swsl2);

                /* position */
                Vector3<scalar> vt = static_cast<scalar>(ksl)*alsl;
                drho[INDRHOX(i,0)] += vt[0];
                drho[INDRHOX(i,1)] += vt[1];
                drho[INDRHOX(i,2)] += vt[2];

                for (si=0; si<R; si++) { /* scale */
                    Vector3<scalar> aisi(rhot[INDRHOP(i,0,si)],rhot[INDRHOP(i,1,si)],rhot[INDRHOP(i,2,si)]);

                    /* scale */
                    vt = (2.0*dot<scalar>(alsl,aisi)*D1Ks(kbasesl,rsl2,swsl2))*ximxl;
                    drho[INDRHOP(i,0,si)] -= vt[0];
                    drho[INDRHOP(i,1,si)] -= vt[1];
                    drho[INDRHOP(i,2,si)] -= vt[2];
                }
            }
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double t;
    double *rhot;
    int L;
    int R;
    int dim;
    double *scales2;
    double *scaleweight2;
    double *drho;

    /* check for proper number of arguments */
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("gradC:nrhs","7 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("gradC:nlhs","One output required.");
    }
    
    /* get the values of the non-array inputs  */
    t = mxGetScalar(prhs[0]);
    L = (double)(mxGetScalar(prhs[2]));
    R = (double)(mxGetScalar(prhs[3]));
    dim = (double)(mxGetScalar(prhs[4]));
    mexAssert(dim == 3);

    int CSP = dim*(1+R);

    /* get the array inputs */
    int n= mxGetM(prhs[1]);
    rhot = mxGetPr(prhs[1]);
    if(CSP*L != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("gradC:nlhs","rhot dimension mismatch (got %d, expected %d).",mxGetM(prhs[1]),CSP*L);
    }
    scales2 = mxGetPr(prhs[5]);
    if(R != mxGetN(prhs[5])) {
        mexErrMsgIdAndTxt("gradC:nlhs","scales2 dimension mismatch (got %d, expected %d).",mxGetN(prhs[5]),R);
    }
    scaleweight2 = mxGetPr(prhs[6]);
    if(R != mxGetN(prhs[6])) {
        mexErrMsgIdAndTxt("gradC:nlhs","scaleweight2 dimension mismatch (got %d, expected %d).",mxGetN(prhs[6]),R);
    }

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    drho = mxGetPr(plhs[0]);

    /* debug */
    /*mexPrintf("%f %f %f %d %d %d %d %f %f %d %d %d\n",t,yt[0],rhot[0],L,R,CSP,CSD,scales2[0],scaleweight2[0],GL,GS,GI);*/

    /* call the computational routine */
    fastPointPathOrder0(drho,t,rhot,L,R,CSP,scales2,scaleweight2);
}
