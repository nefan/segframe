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
void fastPointTransportOrder1(scalar *dgrid, scalar t, scalar *gridt, scalar *Gt, scalar *rhoj0, int L, int R, int CSP, scalar *scales2, scalar *scaleweight2, int Ngrid)
{
    /* dgrid already initialized to zeros */
    const int dim = 3;

#pragma omp parallel for schedule(static) shared(dgrid,t,gridt,Gt,rhoj0,L,R,dim,CSP,scales2,scaleweight2)
    for (int n=0; n<Ngrid; n++) { /* grid point */
        Vector3<scalar> xn(gridt[3*n+0],gridt[3*n+1],gridt[3*n+2]);

        for (int i=0; i<L; i++) { /* particle */
            Vector3<scalar> xi(Gt[INDGX(i,0)],Gt[INDGX(i,1)],Gt[INDGX(i,2)]);
            Vector3<scalar> mui(Gt[INDGMU(i,0)],Gt[INDGMU(i,1)],Gt[INDGMU(i,2)]);
            Matrix3<scalar> DphiiT(
                Vector3<scalar>(Gt[INDGDPHI(i,0,0)],Gt[INDGDPHI(i,1,0)],Gt[INDGDPHI(i,2,0)]),
                Vector3<scalar>(Gt[INDGDPHI(i,0,1)],Gt[INDGDPHI(i,1,1)],Gt[INDGDPHI(i,2,1)]),
                Vector3<scalar>(Gt[INDGDPHI(i,0,2)],Gt[INDGDPHI(i,1,2)],Gt[INDGDPHI(i,2,2)]));
            Matrix3<scalar> Dphii = DphiiT.T();
            Matrix3<scalar> DphiiTinv = DphiiT.inv();

            Vector3<scalar> ximxn = xi-xn;
            scalar v = dot<scalar>(ximxn,ximxn);

            scalar r2 = scales2[0];
            scalar sw2 = scaleweight2[0];
            scalar kbase = exp(-v/r2);

            Vector3<scalar> vt1 = KKs(kbase,sw2)*mui;
            dgrid[3*n+0] += vt1[0];
            dgrid[3*n+1] += vt1[1];
            dgrid[3*n+2] += vt1[2];

            for (int j=0; j<dim; j++) {
                Vector3<scalar> rhoj0j(rhoj0[INDRHOJ(i,0,j)],rhoj0[INDRHOJ(i,1,j)],rhoj0[INDRHOJ(i,2,j)]);
                Vector3<scalar> muij = DphiiTinv*rhoj0j;

                vt1 = (2.0*D1Ks(kbase,r2,sw2)*dot<scalar>(ximxn,DphiiT.row(j)))*muij;
                dgrid[3*n+0] += vt1[0];
                dgrid[3*n+1] += vt1[1];
                dgrid[3*n+2] += vt1[2];
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
    double *Gt;
    double *rhoj0;
    int L;
    int R;
    int dim;
    double *scales2;
    double *scaleweight2;
    double *dgrid;

    /* check for proper number of arguments */
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("gradC:nrhs","9 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("gradC:nlhs","One output required.");
    }
    
    /* get the values of the non-array inputs  */
    t = mxGetScalar(prhs[0]);
    L = (double)(mxGetScalar(prhs[4]));
    R = (double)(mxGetScalar(prhs[5]));
    dim = (double)(mxGetScalar(prhs[6]));
    mexAssert(dim == 3);

    int CSP = dim*(1+dim+1);

    /* get the array inputs */
    gridt = mxGetPr(prhs[1]);
    int n= mxGetM(prhs[1]);
    int Ngrid = n/3;
    Gt = mxGetPr(prhs[2]);
    if(CSP*L != mxGetM(prhs[2])) {
        mexErrMsgIdAndTxt("gradC:nlhs","Gt dimension mismatch (got %d, expected %d).",mxGetM(prhs[2]),CSP*L);
    }
    rhoj0 = mxGetPr(prhs[3]);
    if(dim*dim*L != mxGetM(prhs[3])) {
        mexErrMsgIdAndTxt("gradC:nlhs","rhoj0 dimension mismatch (got %d, expected %d).",mxGetM(prhs[3]),dim*dim*L);
    }
    scales2 = mxGetPr(prhs[7]);
    if(R != mxGetN(prhs[7])) {
        mexErrMsgIdAndTxt("gradC:nlhs","scales2 dimension mismatch (got %d, expected %d).",mxGetN(prhs[7]),R);
    }
    scaleweight2 = mxGetPr(prhs[8]);
    if(R != mxGetN(prhs[8])) {
        mexErrMsgIdAndTxt("gradC:nlhs","scaleweight2 dimension mismatch (got %d, expected %d).",mxGetN(prhs[8]),R);
    }

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    dgrid = mxGetPr(plhs[0]);

    /* debug */
    /*mexPrintf("%f %f %f %d %d %d %d %f %f %d %d %d\n",t,yt[0],rhot[0],L,R,CSP,CSD,scales2[0],scaleweight2[0],GL,GS,GI);*/

    /* call the computational routine */
    fastPointTransportOrder1(dgrid,t,gridt,Gt,rhoj0,L,R,CSP,scales2,scaleweight2,Ngrid);
}
