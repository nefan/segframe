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
#include "order1.hpp"

using namespace sla;
using namespace kd;

/* The computational routine */
void fastPointPathOrder1(scalar *dGt, scalar t, scalar *Gt, scalar *rhoj0, int L, int R, int CSP, scalar *scales2, scalar *scaleweight2)
{
    /* dGt already initialized to zeros */
    const int dim = 3;

#pragma omp parallel for schedule(static) shared(dgrid,t,gridt,Gt,rhoj0,L,R,dim,CSP,scales2,scaleweight2)
    for (int n=0; n<L; n++) { /* grid point */
        Vector3<scalar> xn(Gt[INDGX(n,0)],Gt[INDGX(n,1)],Gt[INDGX(n,2)]);
        Vector3<scalar> mun(Gt[INDGMU(n,0)],Gt[INDGMU(n,1)],Gt[INDGMU(n,2)]);
        Matrix3<scalar> DphinT(
            Vector3<scalar>(Gt[INDGDPHI(n,0,0)],Gt[INDGDPHI(n,1,0)],Gt[INDGDPHI(n,2,0)]),
            Vector3<scalar>(Gt[INDGDPHI(n,0,1)],Gt[INDGDPHI(n,1,1)],Gt[INDGDPHI(n,2,1)]),
            Vector3<scalar>(Gt[INDGDPHI(n,0,2)],Gt[INDGDPHI(n,1,2)],Gt[INDGDPHI(n,2,2)]));
        Matrix3<scalar> DphinTinv = DphinT.inv();

        Vector3<scalar> dphin;
        Matrix3<scalar> dDphinT;
        Vector3<scalar> dmun;

        for (int i=0; i<L; i++) { /* particle */
            Vector3<scalar> xi(Gt[INDGX(i,0)],Gt[INDGX(i,1)],Gt[INDGX(i,2)]);
            Vector3<scalar> mui(Gt[INDGMU(i,0)],Gt[INDGMU(i,1)],Gt[INDGMU(i,2)]);
            Matrix3<scalar> DphiiT(
                Vector3<scalar>(Gt[INDGDPHI(i,0,0)],Gt[INDGDPHI(i,1,0)],Gt[INDGDPHI(i,2,0)]),
                Vector3<scalar>(Gt[INDGDPHI(i,0,1)],Gt[INDGDPHI(i,1,1)],Gt[INDGDPHI(i,2,1)]),
                Vector3<scalar>(Gt[INDGDPHI(i,0,2)],Gt[INDGDPHI(i,1,2)],Gt[INDGDPHI(i,2,2)]));
            //Matrix3<scalar> Dphii = DphiiT.T();
            Matrix3<scalar> DphiiTinv = DphiiT.inv();

            Vector3<scalar> ximxn = xi-xn;
            scalar r2 = scales2[0];
            scalar sw2 = scaleweight2[0];

            // phi
            dphin = dphin + Ks(ximxn,r2,sw2)*mui;

            for (int j=0; j<dim; j++) {
                Vector3<scalar> rhoj0j(rhoj0[INDRHOJ(i,0,j)],rhoj0[INDRHOJ(i,1,j)],rhoj0[INDRHOJ(i,2,j)]);
                Vector3<scalar> muij = DphiiTinv*rhoj0j;

                dphin = dphin + dot<scalar>(N1Ks(ximxn,r2,sw2),DphiiT.row(j))*muij;
            }
            // dPhi
            for (int l=0; l<dim; l++) {
                dDphinT.addRow(l,
                        dot<scalar>(N2Ks(ximxn,r2,sw2),DphinT.row(l))*mui
                        );

                for (int j=0; j<dim; j++) {
                    Vector3<scalar> rhoj0j(rhoj0[INDRHOJ(i,0,j)],rhoj0[INDRHOJ(i,1,j)],rhoj0[INDRHOJ(i,2,j)]);
                    Vector3<scalar> muij = DphiiTinv*rhoj0j;

                    dDphinT.addRow(l,
                            dot<scalar>(D1N2Ks(ximxn,DphiiT.row(j),r2,sw2),DphinT.row(l))*muij
                            );
                }
            }
            // mu
            dmun = dmun - dot<scalar>(mun,mui)*N2Ks(ximxn,r2,sw2);
            for (int j=0; j<dim; j++) {
                Vector3<scalar> rhoj0j(rhoj0[INDRHOJ(i,0,j)],rhoj0[INDRHOJ(i,1,j)],rhoj0[INDRHOJ(i,2,j)]);
                Vector3<scalar> muij = DphiiTinv*rhoj0j;
                Vector3<scalar> rhoj0nj(rhoj0[INDRHOJ(n,0,j)],rhoj0[INDRHOJ(n,1,j)],rhoj0[INDRHOJ(n,2,j)]);
                Vector3<scalar> munj = DphinTinv*rhoj0nj;

                dmun = dmun - dot<scalar>(munj,mui)*D2N2Ks(ximxn,DphinT.row(j),r2,sw2)
                    - dot<scalar>(mun,muij)*D1N2Ks(ximxn,DphiiT.row(j),r2,sw2);

                for (int jm=0; jm<dim; jm++) {
                    Vector3<scalar> rhoj0njm(rhoj0[INDRHOJ(n,0,jm)],rhoj0[INDRHOJ(n,1,jm)],rhoj0[INDRHOJ(n,2,jm)]);
                    Vector3<scalar> munjm = DphinTinv*rhoj0njm;

                    dmun = dmun - dot<scalar>(munjm,muij)*D2D1N2Ks(ximxn,DphiiT.row(j),DphinT.row(jm),r2,sw2);
                }

            }
        }
            
        // store
        dphin.store(&dGt[INDGX(n,0)]);
        dDphinT.store(&dGt[INDGDPHI(n,0,0)]);
        dmun.store(&dGt[INDGMU(n,0)]);
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double t;
    double *Gt;
    double *rhoj0;
    int L;
    int R;
    int dim;
    double *scales2;
    double *scaleweight2;
    double *dGt;

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

    int CSP = dim*(1+dim+1);

    /* get the array inputs */
    int n= mxGetM(prhs[1]);
    Gt = mxGetPr(prhs[1]);
    if(CSP*L != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("gradC:nlhs","Gt dimension mismatch (got %d, expected %d).",mxGetM(prhs[1]),CSP*L);
    }
    rhoj0 = mxGetPr(prhs[2]);
    if(dim*dim*L != mxGetM(prhs[2])) {
        mexErrMsgIdAndTxt("gradC:nlhs","rhoj0 dimension mismatch (got %d, expected %d).",mxGetM(prhs[2]),dim*dim*L);
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
    dGt = mxGetPr(plhs[0]);

    /* debug */
    /*mexPrintf("%f %f %f %d %d %d %d %f %f %d %d %d\n",t,yt[0],rhot[0],L,R,CSP,CSD,scales2[0],scaleweight2[0],GL,GS,GI);*/

    /* call the computational routine */
    fastPointPathOrder1(dGt,t,Gt,rhoj0,L,R,CSP,scales2,scaleweight2);
}
