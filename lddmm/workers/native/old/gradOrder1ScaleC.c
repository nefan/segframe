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

/* energy gradient (helper function) */
void addEgradt(double *dy, double *yt, double *rhot, int L, int R, int dim, int CSP, int CSD, double *scales2, double *scaleweight2, double energyweight)
{
    int CSPL = CSP*L;
    int CSDL = 1;

    int i, l, sl, si;

    if (dim == 2) {
#pragma omp parallel for schedule(static) shared(dy,yt,rhot,L,R,dim,CSP,CSD,CSPL,CSDL,scales2,scaleweight2,energyweight) private(i,l,sl,si)
    for (i=0; i<L; i++) { /* particle */
        double xi[2] = {rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)]};

        /* debug */
        /*mexPrintf("xi %f %f\n",xi[0],xi[1]);*/

        for (l=0; l<L; l++) { /* particle */
            double xl[2] = {rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)]};

            /* debug */
            /*mexPrintf("xl %f %f\n",xl[0],xl[1]);*/

            double ximxl[2]; VECMINUS(ximxl,xi,xl);
            double v = VECDOT2(ximxl,ximxl);

            for (sl=0; sl<R; sl++) { /* scale */
                double alsl[2] = {rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)]};
                double aisl[2] = {rhot[INDRHOP(i,0,sl)],rhot[INDRHOP(i,1,sl)]};
                double rsl2 = scales2[sl];
                double swsl2 = scaleweight2[sl];
                double kbasesl = exp(-v/rsl2);
                double ksl = Ks(kbasesl,swsl2);
                double d1ksl = D1Ks(kbasesl,rsl2,swsl2);
                double d1kslXximxl[2]; VECSCALAR(d1kslXximxl,d1ksl,ximxl);
                double d2ksl = D2Ks(kbasesl,rsl2,swsl2);

                double aislTXalsl; VECDOT(aislTXalsl,aisl,alsl);

                /* dx */
                double vt1[2]; VECSCALAR(vt1,energyweight*4*d1ksl*aislTXalsl,ximxl);
                dy[INDYX(i,0,0,0,0)] += vt1[0];
                dy[INDYX(i,1,0,0,0)] += vt1[1];

                /* da */
                VECSCALAR(vt1,energyweight*2*ksl,alsl);
                dy[INDYP(i,0,sl,0,0,0)] += vt1[0];
                dy[INDYP(i,1,sl,0,0,0)] += vt1[1];
            }
        }
    }
    } else {
#pragma omp parallel for schedule(static) shared(dy,yt,rhot,L,R,dim,CSP,CSD,CSPL,CSDL,scales2,scaleweight2,energyweight) private(i,l,sl,si)
    for (i=0; i<L; i++) { /* particle */
        double xi[3] = {rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)],rhot[INDRHOX(i,2)]};

        for (l=0; l<L; l++) { /* particle */
            double xl[3] = {rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)],rhot[INDRHOX(l,2)]};

            double ximxl[3]; _3VECMINUS(ximxl,xi,xl);
            double v = _3VECDOT2(ximxl,ximxl);

            for (sl=0; sl<R; sl++) { /* scale */
                double alsl[3] = {rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)],rhot[INDRHOP(l,2,sl)]};
                double aisl[3] = {rhot[INDRHOP(i,0,sl)],rhot[INDRHOP(i,1,sl)],rhot[INDRHOP(i,2,sl)]};
                double rsl2 = scales2[sl];
                double swsl2 = scaleweight2[sl];
                double kbasesl = exp(-v/rsl2);
                double ksl = Ks(kbasesl,swsl2);
                double d1ksl = D1Ks(kbasesl,rsl2,swsl2);
                double d1kslXximxl[3]; _3VECSCALAR(d1kslXximxl,d1ksl,ximxl);
                double d2ksl = D2Ks(kbasesl,rsl2,swsl2);

                double aislTXalsl; _3VECDOT(aislTXalsl,aisl,alsl);

                /* dx */
                double vt1[3]; _3VECSCALAR(vt1,energyweight*4*d1ksl*aislTXalsl,ximxl);
                dy[INDYX(i,0,0,0,0)] += vt1[0];
                dy[INDYX(i,1,0,0,0)] += vt1[1];
                dy[INDYX(i,2,0,0,0)] += vt1[2];

                /* da */
                _3VECSCALAR(vt1,energyweight*2*ksl,alsl);
                dy[INDYP(i,0,sl,0,0,0)] += vt1[0];
                dy[INDYP(i,1,sl,0,0,0)] += vt1[1];
                dy[INDYP(i,2,sl,0,0,0)] += vt1[2];
            }
        }
    }
    }
}

/* The computational routine */
void gradOrder1ScaleC(double *dy, double *yt, double *rhot, int L, int R, int dim, int CSP, int CSD, double *scales2, double *scaleweight2, double energyweight)
{
    /* dy already initialized to zeros */
    int CSPL = CSP*L;
    int CSDL = 1;
    double invR = 1.0/R;

    int i, l, sl, si;

    if (dim == 2) {
#pragma omp parallel for schedule(static) shared(dy,yt,rhot,L,R,dim,CSP,CSD,CSPL,CSDL,scales2,scaleweight2) private(i,l,sl,si)
    for (i=0; i<L; i++) { /* particle */
        double xi[2] = {rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)]};
        double dxi[2] = {yt[INDYX(i,0,0,0,0)],yt[INDYX(i,1,0,0,0)]};

        for (l=0; l<L; l++) { /* particle */
            double xl[2] = {rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)]};
            double dxl[2] = {yt[INDYX(l,0,0,0,0)],yt[INDYX(l,1,0,0,0)]};

            double ximxl[2]; VECMINUS(ximxl,xi,xl);
            double v = VECDOT2(ximxl,ximxl);

            for (sl=0; sl<R; sl++) { /* scale */
                double alsl[2] = {rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)]};
                double aisl[2] = {rhot[INDRHOP(i,0,sl)],rhot[INDRHOP(i,1,sl)]};
                double rsl2 = scales2[sl];
                double swsl2 = scaleweight2[sl];
                double kbasesl = exp(-v/rsl2);
                double d1ksl = D1Ks(kbasesl,rsl2,swsl2);
                double d1kslXximxl[2]; VECSCALAR(d1kslXximxl,d1ksl,ximxl);

                double dalsl[2] = {yt[INDYP(l,0,sl,0,0,0)],yt[INDYP(l,1,sl,0,0,0)]};

                /* dx */
                double rt1; VECDOT(rt1,aisl,dxl);
                double rt2; VECDOT(rt2,alsl,dxi);
                double vt1[2]; VECSCALAR(vt1,2*(rt1+rt2),d1kslXximxl);
                dy[INDYX(i,0,0,0,0)] += vt1[0];
                dy[INDYX(i,1,0,0,0)] += vt1[1];

                for (si=0; si<R; si++) { /* scale */
                    double alsi[2] = {rhot[INDRHOP(l,0,si)],rhot[INDRHOP(l,1,si)]};
                    double aisi[2] = {rhot[INDRHOP(i,0,si)],rhot[INDRHOP(i,1,si)]};
                    double daisi[2] = {yt[INDYP(i,0,si,0,0,0)],yt[INDYP(i,1,si,0,0,0)]};
                    double daisl[2] = {yt[INDYP(i,0,sl,0,0,0)],yt[INDYP(i,1,sl,0,0,0)]};

                    double rsi2 = scales2[si];
                    double swsi2 = scaleweight2[si];
                    double kbasesi = exp(-v/rsi2);
                    double ksi = Ks(kbasesi,swsi2);
                    double d1ksi = D1Ks(kbasesi,rsi2,swsi2);
                    double d2ksi = D2Ks(kbasesi,rsi2,swsi2);

                    double aislTXalsi; VECDOT(aislTXalsi,aisl,alsi);
                    double aisiTXalsl; VECDOT(aisiTXalsl,aisi,alsl);

                    /* dx */
                    VECSCALAR(vt1,aislTXalsi,daisl);
                    double vt2[2]; VECSCALAR(vt2,aisiTXalsl,dalsl);
                    double vt3[2]; VECMINUS(vt3,vt1,vt2);
                    double rt; VECDOT(rt,ximxl,vt3);
                    VECSCALAR(vt1,-2*d1ksi,vt3);
                    VECSCALAR(vt2,-4*d2ksi*rt,ximxl);
                    VECADD(vt3,vt1,vt2);
                    dy[INDYX(i,0,0,0,0)] += vt3[0];
                    dy[INDYX(i,1,0,0,0)] += vt3[1];

                    /* da */
                    VECSCALAR(vt1,d1ksl,daisi);
                    VECSCALAR(vt2,d1ksi,dalsl);
                    VECMINUS(vt3,vt1,vt2);
                    VECDOT(rt,ximxl,vt3);
                    VECSCALAR(vt1,-2*rt,alsl);
                    VECSCALAR(vt2,invR*ksi,dxl);
                    VECADD(vt3,vt1,vt2);

                    dy[INDYP(i,0,si,0,0,0)] += vt3[0];
                    dy[INDYP(i,1,si,0,0,0)] += vt3[1];
                }
            }
        }
    }
    } else {
/* #pragma omp parallel for schedule(static) shared(dy,yt,rhot,L,R,dim,CSP,CSD,CSPL,CSDL,scales2,scaleweight2) private(i,l,sl,si)*/
    for (i=0; i<L; i++) { /* particle */
        double xi[3] = {rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)],rhot[INDRHOX(i,2)]};
        double dxi[3] = {yt[INDYX(i,0,0,0,0)],yt[INDYX(i,1,0,0,0)],yt[INDYX(i,2,0,0,0)]};

        for (l=0; l<L; l++) { /* particle */
            double xl[3] = {rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)],rhot[INDRHOX(l,2)]};
            double dxl[3] = {yt[INDYX(l,0,0,0,0)],yt[INDYX(l,1,0,0,0)],yt[INDYX(l,2,0,0,0)]};

            double ximxl[3]; _3VECMINUS(ximxl,xi,xl);
            double v = _3VECDOT2(ximxl,ximxl);

            for (sl=0; sl<R; sl++) { /* scale */
                double alsl[3] = {rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)],rhot[INDRHOP(l,2,sl)]};
                double aisl[3] = {rhot[INDRHOP(i,0,sl)],rhot[INDRHOP(i,1,sl)],rhot[INDRHOP(i,2,sl)]};
                double rsl2 = scales2[sl];
                double swsl2 = scaleweight2[sl];
                double kbasesl = exp(-v/rsl2);
                double d1ksl = D1Ks(kbasesl,rsl2,swsl2);
                double d1kslXximxl[3]; _3VECSCALAR(d1kslXximxl,d1ksl,ximxl);

                double dalsl[3] = {yt[INDYP(l,0,sl,0,0,0)],yt[INDYP(l,1,sl,0,0,0)],yt[INDYP(l,2,sl,0,0,0)]};

                /* dx */
                double rt1; _3VECDOT(rt1,aisl,dxl);
                double rt2; _3VECDOT(rt2,alsl,dxi);
                double vt1[3]; _3VECSCALAR(vt1,2*(rt1+rt2),d1kslXximxl);
                dy[INDYX(i,0,0,0,0)] += vt1[0];
                dy[INDYX(i,1,0,0,0)] += vt1[1];
                dy[INDYX(i,2,0,0,0)] += vt1[2];

                for (si=0; si<R; si++) { /* scale */
                    double alsi[3] = {rhot[INDRHOP(l,0,si)],rhot[INDRHOP(l,1,si)],rhot[INDRHOP(l,2,si)]};
                    double aisi[3] = {rhot[INDRHOP(i,0,si)],rhot[INDRHOP(i,1,si)],rhot[INDRHOP(i,2,si)]};
                    double daisi[3] = {yt[INDYP(i,0,si,0,0,0)],yt[INDYP(i,1,si,0,0,0)],yt[INDYP(i,2,si,0,0,0)]};
                    double daisl[3] = {yt[INDYP(i,0,sl,0,0,0)],yt[INDYP(i,1,sl,0,0,0)],yt[INDYP(i,2,sl,0,0,0)]};

                    double rsi2 = scales2[si];
                    double swsi2 = scaleweight2[si];
                    double kbasesi = exp(-v/rsi2);
                    double ksi = Ks(kbasesi,swsi2);
                    double d1ksi = D1Ks(kbasesi,rsi2,swsi2);
                    double d2ksi = D2Ks(kbasesi,rsi2,swsi2);

                    double aislTXalsi; _3VECDOT(aislTXalsi,aisl,alsi);
                    double aisiTXalsl; _3VECDOT(aisiTXalsl,aisi,alsl);

                    /* dx */
                    _3VECSCALAR(vt1,aislTXalsi,daisl);
                    double vt2[3]; _3VECSCALAR(vt2,aisiTXalsl,dalsl);
                    double vt3[3]; _3VECMINUS(vt3,vt1,vt2);
                    double rt; _3VECDOT(rt,ximxl,vt3);
                    _3VECSCALAR(vt1,-2*d1ksi,vt3);
                    _3VECSCALAR(vt2,-4*d2ksi*rt,ximxl);
                    _3VECADD(vt3,vt1,vt2);
                    dy[INDYX(i,0,0,0,0)] += vt3[0];
                    dy[INDYX(i,1,0,0,0)] += vt3[1];
                    dy[INDYX(i,2,0,0,0)] += vt3[2];

                    /* da */
                    _3VECSCALAR(vt1,d1ksl,daisi);
                    _3VECSCALAR(vt2,d1ksi,dalsl);
                    _3VECMINUS(vt3,vt1,vt2);
                    _3VECDOT(rt,ximxl,vt3);
                    _3VECSCALAR(vt1,-2*rt,alsl);
                    _3VECSCALAR(vt2,invR*ksi,dxl);
                    _3VECADD(vt3,vt1,vt2);

                    dy[INDYP(i,0,si,0,0,0)] += vt3[0];
                    dy[INDYP(i,1,si,0,0,0)] += vt3[1];
                    dy[INDYP(i,2,si,0,0,0)] += vt3[2];
                }
            }
        }
    }
    }

    addEgradt(dy,yt,rhot,L,R,dim,CSP,CSD,scales2,scaleweight2,energyweight);
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
    gradOrder1ScaleC(dy,yt,rhot,L,R,dim,CSP,CSD,scales2,scaleweight2,energyweight);
}
