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
void gradG(double *dy, double t, double *yt, double *rhot, int L, int R, int dim, int CSP, int CSD, double *scales2, double *scaleweight2)
{
    /* dy already initialized to zeros */
    int CSPL = CSP*L;
    int CSDL = CSD*L;

    int i, l, sl, d, ds ,si;

    if (dim == 2) {
#pragma omp parallel for schedule(static) shared(dy,t,yt,rhot,L,R,dim,CSP,CSD,CSPL,CSDL,scales2,scaleweight2) private(i,l,sl,d,ds,si)
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
                double al[2] = {rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)]};
                double rsl2 = scales2[sl];
                double swsl2 = scaleweight2[sl];
                double kbasesl = exp(-v/rsl2);
                double ksl = Ks(kbasesl,swsl2);
                double d1ksl = D1Ks(kbasesl,rsl2,swsl2);
                double d1kslXximxl[2]; VECSCALAR(d1kslXximxl,d1ksl,ximxl);
                double d2ksl = D2Ks(kbasesl,rsl2,swsl2);

                /* debug */
                /*mexPrintf("al %f %f\n",al[0],al[1]);*/

                for (d=0; d<L; d++) { /* direction particle */

                    for (ds=0; ds<R; ds++) { /* direction scale */
                        double dxi[2][2] = {{yt[INDYX(i,0,d,0,ds)],yt[INDYX(i,0,d,1,ds)]},{yt[INDYX(i,1,d,0,ds)],yt[INDYX(i,1,d,1,ds)]}};
                        double dxl[2][2] = {{yt[INDYX(l,0,d,0,ds)],yt[INDYX(l,0,d,1,ds)]},{yt[INDYX(l,1,d,0,ds)],yt[INDYX(l,1,d,1,ds)]}};
                        double dal[2][2] = {{yt[INDYP(l,0,sl,d,0,ds)],yt[INDYP(l,0,sl,d,1,ds)]},{yt[INDYP(l,1,sl,d,0,ds)],yt[INDYP(l,1,sl,d,1,ds)]}};

                        /* debug */
                        /*mexPrintf("dxi %f %f %f %f\n",dxi[0][0],dxi[0][1],dxi[1][0],dxi[1][1]);*/

                        double dximdxl[2][2];
                        MMINUS(dximdxl,dxi,dxl);
                        double vt1[2] = {ximxl[0]*dximdxl[0][0]+ximxl[1]*dximdxl[1][0],ximxl[0]*dximdxl[0][1]+ximxl[1]*dximdxl[1][1]};
                        double ximxlXximxlTXdximdxl[2][2] = {{vt1[0]*ximxl[0],vt1[1]*ximxl[0]},{vt1[0]*ximxl[1],vt1[1]*ximxl[1]}};
                        double TWOd2kslXximxlXximxlTXdximdxl[2][2]; MSCALAR(TWOd2kslXximxlXximxlTXdximdxl,2*d2ksl,ximxlXximxlTXdximdxl);
                        double d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl[2][2]; 
                        MSCALAR(d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl,d1ksl,dximdxl);
                        MADD(d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl,d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl,TWOd2kslXximxlXximxlTXdximdxl);

                        /* dx */
                        double mt1[2][2]; MSCALAR(mt1,ksl,dal);
                        double vt2[2]; VECSCALAR(vt2,2*d1ksl,vt1);
                        double mt2[2][2] = {{vt2[0]*al[0],vt2[1]*al[0]},{vt2[0]*al[1],vt2[1]*al[1]}};
                        MADD(mt1,mt1,mt2);
                        dy[INDYX(i,0,d,0,ds)] += mt1[0][0];
                        dy[INDYX(i,0,d,1,ds)] += mt1[0][1];
                        dy[INDYX(i,1,d,0,ds)] += mt1[1][0];
                        dy[INDYX(i,1,d,1,ds)] += mt1[1][1];

                        for (si=0; si<R; si++) { /* scale */
                            double ai[2] = {rhot[INDRHOP(i,0,si)],rhot[INDRHOP(i,1,si)]};
                            double dai[2][2] = {{yt[INDYP(i,0,si,d,0,ds)],yt[INDYP(i,0,si,d,1,ds)]},{yt[INDYP(i,1,si,d,0,ds)],yt[INDYP(i,1,si,d,1,ds)]}};

                            double aiTXal; VECDOT(aiTXal,ai,al);

                            /* da */
                            MSCALAR(mt1,aiTXal,d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl);
                            vt1[0] = al[0]*dai[0][0]+al[1]*dai[1][0];
                            vt1[1] = al[0]*dai[0][1]+al[1]*dai[1][1];
                            vt2[0] = ai[0]*dal[0][0]+ai[1]*dal[1][0];
                            vt2[1] = ai[0]*dal[0][1]+ai[1]*dal[1][1];
                            VECADD(vt1,vt1,vt2);
                            mt2[0][0] = vt1[0]*d1kslXximxl[0]; mt2[1][0] = vt1[0]*d1kslXximxl[1];
                            mt2[0][1] = vt1[1]*d1kslXximxl[0]; mt2[1][1] = vt1[1]*d1kslXximxl[1];

                            MADD(mt1,mt1,mt2);
                            MSCALAR(mt1,-2,mt1);

                            dy[INDYP(i,0,si,d,0,ds)] += mt1[0][0];
                            dy[INDYP(i,0,si,d,1,ds)] += mt1[0][1];
                            dy[INDYP(i,1,si,d,0,ds)] += mt1[1][0];
                            dy[INDYP(i,1,si,d,1,ds)] += mt1[1][1];

                            /* debug */
                            /*mexPrintf("ai %f %f\n",ai[0],ai[1]);*/
                        }
                    }
                }
            }
        }
    }
    } else {
#pragma omp parallel for schedule(static) shared(dy,t,yt,rhot,L,R,dim,CSP,CSD,CSPL,CSDL,scales2,scaleweight2) private(i,l,sl,d,ds,si)
    for (i=0; i<L; i++) { /* particle */
        /*double xi[2] = {rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)]};*/
        double xi[3] = {rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)],rhot[INDRHOX(i,2)]};

        for (l=0; l<L; l++) { /* particle */
            /*double xl[2] = {rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)]};*/
            double xl[3] = {rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)],rhot[INDRHOX(l,2)]};

            /*double ximxl[2]; VECMINUS(ximxl,xi,xl);
            double v = VECDOT2(ximxl,ximxl);*/
            double ximxl[3]; _3VECMINUS(ximxl,xi,xl);
            double v = _3VECDOT2(ximxl,ximxl);

            for (sl=0; sl<R; sl++) { /* scale */
                /*double al[2] = {rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)]};*/
                double al[3] = {rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)],rhot[INDRHOP(l,2,sl)]};
                double rsl2 = scales2[sl];
                double swsl2 = scaleweight2[sl];
                double kbasesl = exp(-v/rsl2);
                double ksl = Ks(kbasesl,swsl2);
                double d1ksl = D1Ks(kbasesl,rsl2,swsl2);
                /*double d1kslXximxl[2]; VECSCALAR(d1kslXximxl,d1ksl,ximxl);*/
                double d1kslXximxl[3]; _3VECSCALAR(d1kslXximxl,d1ksl,ximxl);
                double d2ksl = D2Ks(kbasesl,rsl2,swsl2);

                for (d=0; d<L; d++) { /* direction particle */

                    for (ds=0; ds<R; ds++) { /* direction scale */
                        /*double dxi[2][2] = {{yt[INDYX(i,0,d,0,ds)],yt[INDYX(i,0,d,1,ds)]},{yt[INDYX(i,1,d,0,ds)],yt[INDYX(i,1,d,1,ds)]}};
                        double dxl[2][2] = {{yt[INDYX(l,0,d,0,ds)],yt[INDYX(l,0,d,1,ds)]},{yt[INDYX(l,1,d,0,ds)],yt[INDYX(l,1,d,1,ds)]}};
                        double dal[2][2] = {{yt[INDYP(l,0,sl,d,0,ds)],yt[INDYP(l,0,sl,d,1,ds)]},{yt[INDYP(l,1,sl,d,0,ds)],yt[INDYP(l,1,sl,d,1,ds)]}};*/
                        double dxi[3][3] = {{yt[INDYX(i,0,d,0,ds)],yt[INDYX(i,0,d,1,ds)],yt[INDYX(i,0,d,2,ds)]},{yt[INDYX(i,1,d,0,ds)],yt[INDYX(i,1,d,1,ds)],yt[INDYX(i,1,d,2,ds)]},{yt[INDYX(i,2,d,0,ds)],yt[INDYX(i,2,d,1,ds)],yt[INDYX(i,2,d,2,ds)]}};
                        double dxl[3][3] = {{yt[INDYX(l,0,d,0,ds)],yt[INDYX(l,0,d,1,ds)],yt[INDYX(l,0,d,2,ds)]},{yt[INDYX(l,1,d,0,ds)],yt[INDYX(l,1,d,1,ds)],yt[INDYX(l,1,d,2,ds)]},{yt[INDYX(l,2,d,0,ds)],yt[INDYX(l,2,d,1,ds)],yt[INDYX(l,2,d,2,ds)]}};
                        double dal[3][3] = {{yt[INDYP(l,0,sl,d,0,ds)],yt[INDYP(l,0,sl,d,1,ds)],yt[INDYP(l,0,sl,d,2,ds)]},{yt[INDYP(l,1,sl,d,0,ds)],yt[INDYP(l,1,sl,d,1,ds)],yt[INDYP(l,1,sl,d,2,ds)]},{yt[INDYP(l,2,sl,d,0,ds)],yt[INDYP(l,2,sl,d,1,ds)],yt[INDYP(l,2,sl,d,2,ds)]}};

                        /*double dximdxl[2][2];
                        MMINUS(dximdxl,dxi,dxl);
                        double vt1[2] = {ximxl[0]*dximdxl[0][0]+ximxl[1]*dximdxl[1][0],ximxl[0]*dximdxl[0][1]+ximxl[1]*dximdxl[1][1]};
                        double ximxlXximxlTXdximdxl[2][2] = {{vt1[0]*ximxl[0],vt1[1]*ximxl[0]},{vt1[0]*ximxl[1],vt1[1]*ximxl[1]}};
                        double TWOd2kslXximxlXximxlTXdximdxl[2][2]; MSCALAR(TWOd2kslXximxlXximxlTXdximdxl,2*d2ksl,ximxlXximxlTXdximdxl);
                        double d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl[2][2]; 
                        MSCALAR(d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl,d1ksl,dximdxl);
                        MADD(d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl,d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl,TWOd2kslXximxlXximxlTXdximdxl);*/
                        double dximdxl[3][3];
                        _3MMINUS(dximdxl,dxi,dxl);
                        double vt1[3] = {ximxl[0]*dximdxl[0][0]+ximxl[1]*dximdxl[1][0]+ximxl[2]*dximdxl[2][0],ximxl[0]*dximdxl[0][1]+ximxl[1]*dximdxl[1][1]+ximxl[2]*dximdxl[2][1],ximxl[0]*dximdxl[0][2]+ximxl[1]*dximdxl[1][2]+ximxl[2]*dximdxl[2][2]};
                        double ximxlXximxlTXdximdxl[3][3] = {{vt1[0]*ximxl[0],vt1[1]*ximxl[0],vt1[2]*ximxl[0]},{vt1[0]*ximxl[1],vt1[1]*ximxl[1],vt1[2]*ximxl[1]},{vt1[0]*ximxl[2],vt1[1]*ximxl[2],vt1[2]*ximxl[2]}};
                        double TWOd2kslXximxlXximxlTXdximdxl[3][3]; _3MSCALAR(TWOd2kslXximxlXximxlTXdximdxl,2*d2ksl,ximxlXximxlTXdximdxl);
                        double d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl[3][3]; 
                        _3MSCALAR(d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl,d1ksl,dximdxl);
                        _3MADD(d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl,d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl,TWOd2kslXximxlXximxlTXdximdxl);

                        /* dx */
                        /*double mt1[2][2]; MSCALAR(mt1,ksl,dal);
                        double vt2[2]; VECSCALAR(vt2,2*d1ksl,vt1);
                        double mt2[2][2] = {{vt2[0]*al[0],vt2[1]*al[0]},{vt2[0]*al[1],vt2[1]*al[1]}};
                        MADD(mt1,mt1,mt2);*/
                        double mt1[3][3]; _3MSCALAR(mt1,ksl,dal);
                        double vt2[3]; _3VECSCALAR(vt2,2*d1ksl,vt1);
                        double mt2[3][3] = {{vt2[0]*al[0],vt2[1]*al[0],vt2[2]*al[0]},{vt2[0]*al[1],vt2[1]*al[1],vt2[2]*al[1]},{vt2[0]*al[2],vt2[1]*al[2],vt2[2]*al[2]}};
                        _3MADD(mt1,mt1,mt2);
                        /*dy[INDYX(i,0,d,0,ds)] += mt1[0][0];
                        dy[INDYX(i,0,d,1,ds)] += mt1[0][1];
                        dy[INDYX(i,1,d,0,ds)] += mt1[1][0];
                        dy[INDYX(i,1,d,1,ds)] += mt1[1][1];*/
                        dy[INDYX(i,0,d,0,ds)] += mt1[0][0]; dy[INDYX(i,0,d,1,ds)] += mt1[0][1]; dy[INDYX(i,0,d,2,ds)] += mt1[0][2];
                        dy[INDYX(i,1,d,0,ds)] += mt1[1][0]; dy[INDYX(i,1,d,1,ds)] += mt1[1][1]; dy[INDYX(i,1,d,2,ds)] += mt1[1][2];
                        dy[INDYX(i,2,d,0,ds)] += mt1[2][0]; dy[INDYX(i,2,d,1,ds)] += mt1[2][1]; dy[INDYX(i,2,d,2,ds)] += mt1[2][2];

                        for (si=0; si<R; si++) { /* scale */
                            /*double ai[2] = {rhot[INDRHOP(i,0,si)],rhot[INDRHOP(i,1,si)]};
                            double dai[2][2] = {{yt[INDYP(i,0,si,d,0,ds)],yt[INDYP(i,0,si,d,1,ds)]},{yt[INDYP(i,1,si,d,0,ds)],yt[INDYP(i,1,si,d,1,ds)]}};*/
                            double ai[3] = {rhot[INDRHOP(i,0,si)],rhot[INDRHOP(i,1,si)],rhot[INDRHOP(i,2,si)]};
                            double dai[3][3] = {{yt[INDYP(i,0,si,d,0,ds)],yt[INDYP(i,0,si,d,1,ds)],yt[INDYP(i,0,si,d,2,ds)]},{yt[INDYP(i,1,si,d,0,ds)],yt[INDYP(i,1,si,d,1,ds)],yt[INDYP(i,1,si,d,2,ds)]},{yt[INDYP(i,2,si,d,0,ds)],yt[INDYP(i,2,si,d,1,ds)],yt[INDYP(i,2,si,d,2,ds)]}};

                            /*double aiTXal; VECDOT(aiTXal,ai,al);*/
                            double aiTXal; _3VECDOT(aiTXal,ai,al);

                            /* da */
                            /*MSCALAR(mt1,aiTXal,d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl);
                            vt1[0] = al[0]*dai[0][0]+al[1]*dai[1][0];
                            vt1[1] = al[0]*dai[0][1]+al[1]*dai[1][1];
                            vt2[0] = ai[0]*dal[0][0]+ai[1]*dal[1][0];
                            vt2[1] = ai[0]*dal[0][1]+ai[1]*dal[1][1];
                            VECADD(vt1,vt1,vt2);
                            mt2[0][0] = vt1[0]*d1kslXximxl[0]; mt2[1][0] = vt1[0]*d1kslXximxl[1];
                            mt2[0][1] = vt1[1]*d1kslXximxl[0]; mt2[1][1] = vt1[1]*d1kslXximxl[1];

                            MADD(mt1,mt1,mt2);
                            MSCALAR(mt1,-2,mt1);

                            dy[INDYP(i,0,si,d,0,ds)] += mt1[0][0];
                            dy[INDYP(i,0,si,d,1,ds)] += mt1[0][1];
                            dy[INDYP(i,1,si,d,0,ds)] += mt1[1][0];
                            dy[INDYP(i,1,si,d,1,ds)] += mt1[1][1];*/
                            _3MSCALAR(mt1,aiTXal,d1kslXdximdxlP2d2kslXximxlXximxlTXdximdxl);
                            vt1[0] = al[0]*dai[0][0]+al[1]*dai[1][0]+al[2]*dai[2][0];
                            vt1[1] = al[0]*dai[0][1]+al[1]*dai[1][1]+al[2]*dai[2][1];
                            vt1[2] = al[0]*dai[0][2]+al[1]*dai[1][2]+al[2]*dai[2][2];
                            vt2[0] = ai[0]*dal[0][0]+ai[1]*dal[1][0]+ai[2]*dal[2][0];
                            vt2[1] = ai[0]*dal[0][1]+ai[1]*dal[1][1]+ai[2]*dal[2][1];
                            vt2[2] = ai[0]*dal[0][2]+ai[1]*dal[1][2]+ai[2]*dal[2][2];
                            _3VECADD(vt1,vt1,vt2);
                            mt2[0][0] = vt1[0]*d1kslXximxl[0]; mt2[1][0] = vt1[0]*d1kslXximxl[1]; mt2[2][0] = vt1[0]*d1kslXximxl[2]; 
                            mt2[0][1] = vt1[1]*d1kslXximxl[0]; mt2[1][1] = vt1[1]*d1kslXximxl[1]; mt2[2][1] = vt1[1]*d1kslXximxl[2];
                            mt2[0][2] = vt1[2]*d1kslXximxl[0]; mt2[1][2] = vt1[2]*d1kslXximxl[1]; mt2[2][2] = vt1[2]*d1kslXximxl[2];

                            _3MADD(mt1,mt1,mt2);
                            _3MSCALAR(mt1,-2,mt1);

                            dy[INDYP(i,0,si,d,0,ds)] += mt1[0][0]; dy[INDYP(i,0,si,d,1,ds)] += mt1[0][1]; dy[INDYP(i,0,si,d,2,ds)] += mt1[0][2];
                            dy[INDYP(i,1,si,d,0,ds)] += mt1[1][0]; dy[INDYP(i,1,si,d,1,ds)] += mt1[1][1]; dy[INDYP(i,1,si,d,2,ds)] += mt1[1][2];
                            dy[INDYP(i,2,si,d,0,ds)] += mt1[2][0]; dy[INDYP(i,2,si,d,1,ds)] += mt1[2][1]; dy[INDYP(i,2,si,d,2,ds)] += mt1[2][2];
                        }
                    }
                }
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
    double *yt;
    double *rhot;
    int L;
    int R;
    int dim;
    int CSP;
    int CSD;
    double *scales2;
    double *scaleweight2;
    double *dy;

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
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    dy = mxGetPr(plhs[0]);

    /* debug */
    /*mexPrintf("%f %f %f %d %d %d %d %f %f\n",t,yt[0],rhot[0],L,R,CSP,CSD,scales2[0],scaleweight2[0]);*/

    /* call the computational routine */
    gradG(dy,t,yt,rhot,L,R,dim,CSP,CSD,scales2,scaleweight2);
}
