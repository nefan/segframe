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

/* The computational routine */
void gradEH(double *_w, double t, double *yt, double *rhot, int L, int R, int dim, int CSP, int CSD, double *scales2, double *scaleweight2, int GL, int GSL, int GLi, int GI, int GSI, int GIi)
{
    double w = 0;
    int CSPL = CSP*L;
    int CSDL = CSD*L;

    int i, l, s;

    if (dim == 2) {
    for (i=0; i<L; i++) { /* particle */
        double xi[2] = {rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)]};
        double dxis[2] = {yt[INDYX(i,0,GL,GLi,GSL)],yt[INDYX(i,1,GL,GLi,GSL)]};
        double dxir[2] = {yt[INDYX(i,0,GI,GIi,GSI)],yt[INDYX(i,1,GI,GIi,GSI)]};

        /* debug */
        /*mexPrintf("xi %f %f\n",xi[0],xi[1]);
        mexPrintf("dxi %f %f\n",dxi[0],dxi[1]);*/
        
        for (l=0; l<L; l++) { /* particle */
            double xl[2] = {rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)]};
            double dxls[2] = {yt[INDYX(l,0,GL,GLi,GSL)],yt[INDYX(l,1,GL,GLi,GSL)]};
            double dxlr[2] = {yt[INDYX(l,0,GI,GIi,GSI)],yt[INDYX(l,1,GI,GIi,GSI)]};

            /* debug */
            /*mexPrintf("xl %f %f\n",xl[0],xl[1]);
            mexPrintf("dxl %f %f\n",dxl[0],dxl[1]);*/

            double ximxl[2]; VECMINUS(ximxl,xi,xl);
            double v = VECDOT2(ximxl,ximxl);
            double dxismdxls[2];
            VECMINUS(dxismdxls,dxis,dxls);
            double dxirmdxlr[2];
            VECMINUS(dxirmdxlr,dxir,dxlr);

            for (s=0; s<R; s++) { /* scale */
                double ai[2] = {rhot[INDRHOP(i,0,s)],rhot[INDRHOP(i,1,s)]};
                double al[2] = {rhot[INDRHOP(l,0,s)],rhot[INDRHOP(l,1,s)]};

                double dals[2] = {yt[INDYP(l,0,s,GL,GLi,GSL)],yt[INDYP(l,1,s,GL,GLi,GSL)]};
                double dair[2] = {yt[INDYP(i,0,s,GI,GIi,GSI)],yt[INDYP(i,1,s,GI,GIi,GSI)]};
                double dalr[2] = {yt[INDYP(l,0,s,GI,GIi,GSI)],yt[INDYP(l,1,s,GI,GIi,GSI)]};

                /* debug */
                /*mexPrintf("ai %f %f\n",ai[0],ai[1]);
                mexPrintf("al %f %f\n",al[0],al[1]);
                mexPrintf("dai %f %f\n",dai[0],dai[1]);
                mexPrintf("dal %f %f\n",dal[0],dal[1]);*/

                double r2 = scales2[s];
                double sw2 = scaleweight2[s];
                double kbase = exp(-v/r2);
                double ks = Ks(kbase,sw2);
                double d1k = D1Ks(kbase,r2,sw2);
                double d1kXximxl[2]; VECSCALAR(d1kXximxl,d1k,ximxl);
                double d2k = D2Ks(kbase,r2,sw2);

                double aidal = VECDOT2(ai,al);

                w += 4*(VECDOT2(d1kXximxl,dxirmdxlr)*VECDOT2(ai,dals)+VECDOT2(d1kXximxl,dxismdxls)*VECDOT2(ai,dalr));
                w += 2*ks*VECDOT2(dair,dals);
                w += 4*d2k*VECDOT2(ximxl,dxirmdxlr)*aidal*VECDOT2(ximxl,dxismdxls);
                w += 2*d1k*aidal*VECDOT2(dxirmdxlr,dxismdxls);

                /* debug */
                /*mexPrintf("w %f %f %f %f %f\n",w,ks*VECDOT2(al,dai),VECDOT2(ai,vt1),ks,d1k);*/
            }
        }
    }
    } else {
#pragma omp parallel for schedule(static) shared(w,t,yt,rhot,L,R,dim,CSP,CSD,CSPL,scales2,scaleweight2,GL,GS,GI) private(i,l,s)
#ifdef FALSE
    for (i=0; i<L; i++) { /* particle */
        /*double xi[2] = {rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)]};
        double dxi[2] = {yt[INDYX(i,0,GL,GI,GS)],yt[INDYX(i,1,GL,GI,GS)]};*/
        double xi[3] = {rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)],rhot[INDRHOX(i,2)]};
        double dxi[3] = {yt[INDYX(i,0,GL,GI,GS)],yt[INDYX(i,1,GL,GI,GS)],yt[INDYX(i,2,GL,GI,GS)]};

        for (l=0; l<L; l++) { /* particle */
            /*double xl[2] = {rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)]};
            double dxl[2] = {yt[INDYX(l,0,GL,GI,GS)],yt[INDYX(l,1,GL,GI,GS)]};*/
            double xl[3] = {rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)],rhot[INDRHOX(l,2)]};
            double dxl[3] = {yt[INDYX(l,0,GL,GI,GS)],yt[INDYX(l,1,GL,GI,GS)],yt[INDYX(l,2,GL,GI,GS)]};

            /*double ximxl[2]; VECMINUS(ximxl,xi,xl);
            double v = VECDOT2(ximxl,ximxl);
            double dximdxl[2];
            VECMINUS(dximdxl,dxi,dxl);*/
            double ximxl[3]; _3VECMINUS(ximxl,xi,xl);
            double v = _3VECDOT2(ximxl,ximxl);
            double dximdxl[3];
            _3VECMINUS(dximdxl,dxi,dxl);

            for (s=0; s<R; s++) { /* scale */
                /*double ai[2] = {rhot[INDRHOP(i,0,s)],rhot[INDRHOP(i,1,s)]};
                double al[2] = {rhot[INDRHOP(l,0,s)],rhot[INDRHOP(l,1,s)]};*/
                double ai[3] = {rhot[INDRHOP(i,0,s)],rhot[INDRHOP(i,1,s)],rhot[INDRHOP(i,2,s)]};
                double al[3] = {rhot[INDRHOP(l,0,s)],rhot[INDRHOP(l,1,s)],rhot[INDRHOP(l,2,s)]};

                /*double dai[2] = {yt[INDYP(i,0,s,GL,GI,GS)],yt[INDYP(i,1,s,GL,GI,GS)]};
                double dal[2] = {yt[INDYP(l,0,s,GL,GI,GS)],yt[INDYP(l,1,s,GL,GI,GS)]};*/
                double dai[3] = {yt[INDYP(i,0,s,GL,GI,GS)],yt[INDYP(i,1,s,GL,GI,GS)],yt[INDYP(i,2,s,GL,GI,GS)]};
                double dal[3] = {yt[INDYP(l,0,s,GL,GI,GS)],yt[INDYP(l,1,s,GL,GI,GS)],yt[INDYP(l,2,s,GL,GI,GS)]};

                double r2 = scales2[s];
                double sw2 = scaleweight2[s];
                double kbase = exp(-v/r2);
                double ks = Ks(kbase,sw2);
                double d1k = D1Ks(kbase,r2,sw2);
                /*double d1kXximxl[2]; VECSCALAR(d1kXximxl,d1k,ximxl);*/
                double d1kXximxl[3]; _3VECSCALAR(d1kXximxl,d1k,ximxl);

                /*w += ks*VECDOT2(al,dai);
                double vt1[2]; VECSCALAR(vt1,2*VECDOT2(d1kXximxl,dximdxl),al);
                double vt2[2]; VECSCALAR(vt2,ks,dal);
                VECADD(vt1,vt1,vt2);
                w += VECDOT2(ai,vt1);*/
                w += ks*_3VECDOT2(al,dai);
                double vt1[3]; _3VECSCALAR(vt1,2*_3VECDOT2(d1kXximxl,dximdxl),al);
                double vt2[3]; _3VECSCALAR(vt2,ks,dal);
                _3VECADD(vt1,vt1,vt2);
                w += _3VECDOT2(ai,vt1);
            }
        }
    }
#endif
    }

    /* save result */
    _w[0] = w;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double t;
    double *yt;
    double *rhot;
    int L;
    int R;
    int dim;
    double *scales2;
    double *scaleweight2;
    int GL;
    int GSL;
    int GLi;
    int GI;
    int GSI;
    int GIi;
    int CSP;
    int CSD;
    double *w;

    /* check for proper number of arguments */
    if(nrhs!=14) {
        mexErrMsgIdAndTxt("gradC:nrhs","14 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("gradC:nlhs","One output required.");
    }
    
    /* get the values of the non-array inputs  */
    t = mxGetScalar(prhs[0]);
    L = (double)(mxGetScalar(prhs[3]));
    R = (double)(mxGetScalar(prhs[4]));
    dim = (double)(mxGetScalar(prhs[5]));
    GL = (double)(mxGetScalar(prhs[8]))-1;
    GSL = (double)(mxGetScalar(prhs[9]))-1;
    GLi = (double)(mxGetScalar(prhs[10]))-1;
    GI = (double)(mxGetScalar(prhs[11]))-1;
    GSI = (double)(mxGetScalar(prhs[12]))-1;
    GIi = (double)(mxGetScalar(prhs[13]))-1;

    CSP = dim*(1+R);
    CSD = dim*R;

    /* get the array inputs */
    yt = mxGetPr(prhs[1]);
    int n= mxGetM(prhs[1]);
    if(n != CSP*L*CSD*L) {
        mexErrMsgIdAndTxt("gradC:nlhs","yt dimension mismatch (got %d, expected %d).",n,CSP*L*CSD*L);
    }
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
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    w = mxGetPr(plhs[0]);

    /* debug */
    /*mexPrintf("%f %f %f %d %d %d %d %f %f %d %d %d\n",t,yt[0],rhot[0],L,R,CSP,CSD,scales[0],scaleweight[0],GL,GS,GI);*/

    /* call the computational routine */
    gradEH(w,t,yt,rhot,L,R,dim,CSP,CSD,scales2,scaleweight2,GL,GSL,GLi,GI,GSI,GIi);
}
