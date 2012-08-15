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

/* energy gradient (helper function) */
void addEgradt(scalar *dy, scalar *yt, scalar *rhot, int L, int R, int CSP, int CSD, scalar *scales2, scalar *scaleweight2, scalar energyweight)
{
    int CSPL = CSP*L;
    int CSDL = 1;
    const int dim = 3;

#pragma omp parallel for schedule(static) shared(dy,yt,rhot,L,R,dim,CSP,CSD,CSPL,CSDL,scales2,scaleweight2,energyweight)
    for (int i=0; i<L; i++) { /* particle */
        Vector3<scalar> xi(rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)],rhot[INDRHOX(i,2)]);

        for (int l=0; l<L; l++) { /* particle */
            Vector3<scalar> xl(rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)],rhot[INDRHOX(l,2)]);

            Vector3<scalar> ximxl = xi-xl;
            scalar v = dot<scalar>(ximxl,ximxl);

            for (int sl=0; sl<R; sl++) { /* scale */
                Vector3<scalar> alsl(rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)],rhot[INDRHOP(l,2,sl)]);
                Vector3<scalar> aisl(rhot[INDRHOP(i,0,sl)],rhot[INDRHOP(i,1,sl)],rhot[INDRHOP(i,2,sl)]);
                scalar rsl2 = scales2[sl];
                scalar swsl2 = scaleweight2[sl];
                scalar kbasesl = exp(-v/rsl2);
                scalar ksl = KKs(kbasesl,swsl2);
                scalar d1ksl = D1Ks(kbasesl,rsl2,swsl2);

                /* dx */
                Vector3<scalar> vt1 = (4.0*energyweight*d1ksl*dot<scalar>(aisl,alsl))*ximxl;
                dy[INDYX(i,0,0,0,0)] += vt1[0];
                dy[INDYX(i,1,0,0,0)] += vt1[1];
                dy[INDYX(i,2,0,0,0)] += vt1[2];

                /* da */
                vt1 = (2.0*energyweight*ksl)*alsl;
                dy[INDYP(i,0,sl,0,0,0)] += vt1[0];
                dy[INDYP(i,1,sl,0,0,0)] += vt1[1];
                dy[INDYP(i,2,sl,0,0,0)] += vt1[2];
            }
        }
    }
}

/* The computational routine */
void fastPointGradTransportOrder0(scalar *dy, scalar *yt, scalar *rhot, int L, int R, int CSP, int CSD, scalar *scales2, scalar *scaleweight2, scalar energyweight)
{
    /* dy already initialized to zeros */
    int CSPL = CSP*L;
    int CSDL = 1;
    const int dim = 3;

    scalar invR = 1.0/R;

#pragma omp parallel for schedule(static) shared(dy,yt,rhot,L,R,dim,CSP,CSD,CSPL,CSDL,scales2,scaleweight2)
    for (int i=0; i<L; i++) { /* particle */
        Vector3<scalar> xi(rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)],rhot[INDRHOX(i,2)]);
        Vector3<scalar> dxi(yt[INDYX(i,0,0,0,0)],yt[INDYX(i,1,0,0,0)],yt[INDYX(i,2,0,0,0)]);

        for (int l=0; l<L; l++) { /* particle */
            Vector3<scalar> xl(rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)],rhot[INDRHOX(l,2)]);
            Vector3<scalar> dxl(yt[INDYX(l,0,0,0,0)],yt[INDYX(l,1,0,0,0)],yt[INDYX(l,2,0,0,0)]);

            Vector3<scalar> ximxl = xi-xl;
            scalar v = dot<scalar>(ximxl,ximxl);

            for (int sl=0; sl<R; sl++) { /* scale */
                Vector3<scalar> alsl(rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)],rhot[INDRHOP(l,2,sl)]);
                Vector3<scalar> aisl(rhot[INDRHOP(i,0,sl)],rhot[INDRHOP(i,1,sl)],rhot[INDRHOP(i,2,sl)]);
                Vector3<scalar> dalsl(yt[INDYP(l,0,sl,0,0,0)],yt[INDYP(l,1,sl,0,0,0)],yt[INDYP(l,2,sl,0,0,0)]);
                scalar rsl2 = scales2[sl];
                scalar swsl2 = scaleweight2[sl];
                scalar kbasesl = exp(-v/rsl2);
                scalar d1ksl = D1Ks(kbasesl,rsl2,swsl2);

                /* dx */
                Vector3<scalar> vt1 = (2.0*(dot<scalar>(aisl,dxl)+dot<scalar>(alsl,dxi))*d1ksl)*ximxl;
                dy[INDYX(i,0,0,0,0)] += vt1[0];
                dy[INDYX(i,1,0,0,0)] += vt1[1];
                dy[INDYX(i,2,0,0,0)] += vt1[2];

                for (int si=0; si<R; si++) { /* scale */
                    Vector3<scalar> alsi(rhot[INDRHOP(l,0,si)],rhot[INDRHOP(l,1,si)],rhot[INDRHOP(l,2,si)]);
                    Vector3<scalar> aisi(rhot[INDRHOP(i,0,si)],rhot[INDRHOP(i,1,si)],rhot[INDRHOP(i,2,si)]);
                    Vector3<scalar> daisi(yt[INDYP(i,0,si,0,0,0)],yt[INDYP(i,1,si,0,0,0)],yt[INDYP(i,2,si,0,0,0)]);
                    Vector3<scalar> daisl(yt[INDYP(i,0,sl,0,0,0)],yt[INDYP(i,1,sl,0,0,0)],yt[INDYP(i,2,sl,0,0,0)]);

                    scalar rsi2 = scales2[si];
                    scalar swsi2 = scaleweight2[si];
                    scalar kbasesi = exp(-v/rsi2);
                    scalar ksi = KKs(kbasesi,swsi2);
                    scalar d1ksi = D1Ks(kbasesi,rsi2,swsi2);
                    scalar d2ksi = D2Ks(kbasesi,rsi2,swsi2);

                    /* dx */
                    Vector3<scalar> vt2 = dot<scalar>(aisl,alsi)*daisl-dot<scalar>(aisi,alsl)*dalsl;
                    vt1 = -2.0*d1ksi*vt2-4.0*d2ksi*dot<scalar>(ximxl,vt2)*ximxl;
                    dy[INDYX(i,0,0,0,0)] += vt1[0];
                    dy[INDYX(i,1,0,0,0)] += vt1[1];
                    dy[INDYX(i,2,0,0,0)] += vt1[2];

                    /* da */
                    vt1 = -2.0*dot<scalar>(d1ksl*daisi-d1ksi*dalsl,ximxl)*alsl+invR*ksi*dxl;

                    dy[INDYP(i,0,si,0,0,0)] += vt1[0];
                    dy[INDYP(i,1,si,0,0,0)] += vt1[1];
                    dy[INDYP(i,2,si,0,0,0)] += vt1[2];
                }
            }
        }
    }

    addEgradt(dy,yt,rhot,L,R,CSP,CSD,scales2,scaleweight2,energyweight);
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *yt;
    double *rhot;
    int L;
    int R;
    int dim;
    int CSP;
    int CSD;
    double *scales2;
    double *scaleweight2;
    double energyweight;
    double *dy;

    /* check for proper number of arguments */
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("gradC:nrhs","8 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("gradC:nlhs","One output required.");
    }
    
    /* get the values of the non-array inputs  */
    L = (double)(mxGetScalar(prhs[2]));
    R = (double)(mxGetScalar(prhs[3]));
    dim = (double)(mxGetScalar(prhs[4]));
    mexAssert(dim == 3);
    energyweight = (double)(mxGetScalar(prhs[7]));

    CSP = dim*(1+R);
    CSD = dim*R;

    /* get the array inputs */
    yt = mxGetPr(prhs[0]);
    int n= mxGetM(prhs[0]);
    if(n != CSP*L) {
        mexErrMsgIdAndTxt("gradC:nlhs","yt dimension mismatch (got %d, expected %d).",n,CSP*L);
    }
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
    dy = mxGetPr(plhs[0]);

    /* debug */
    /*mexPrintf("%f %f %f %d %d %d %d %f %f\n",t,yt[0],rhot[0],L,R,CSP,CSD,scales2[0],scaleweight2[0]);*/

    /* call the computational routine */
    fastPointGradTransportOrder0(dy,yt,rhot,L,R,CSP,CSD,scales2,scaleweight2,energyweight);
}
