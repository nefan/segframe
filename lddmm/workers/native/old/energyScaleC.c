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
double energyScaleC(double *Et, double *rhot, int L, int R, int dim, int CSP, double *scales2, double *scaleweight2)
{
    /* Et already initialized to zeros */
    int CSPL = CSP*L;

    int i, l, sl;
    double Ett = 0;

    if (dim == 2) {
#pragma omp parallel for schedule(static) shared(rhot,L,R,dim,CSP,CSPL,scales2,scaleweight2) private(i,l,sl) reduction(+:Ett)
    for (i=0; i<L; i++) { /* particle */
        double xi[2] = {rhot[INDRHOX(i,0)],rhot[INDRHOX(i,1)]};

        for (l=0; l<L; l++) { /* particle */
            double xl[2] = {rhot[INDRHOX(l,0)],rhot[INDRHOX(l,1)]};

            double ximxl[2]; VECMINUS(ximxl,xi,xl);
            double v = VECDOT2(ximxl,ximxl);

            for (sl=0; sl<R; sl++) { /* scale */
                double alsl[2] = {rhot[INDRHOP(l,0,sl)],rhot[INDRHOP(l,1,sl)]};
                double aisl[2] = {rhot[INDRHOP(i,0,sl)],rhot[INDRHOP(i,1,sl)]};
                double rsl2 = scales2[sl];
                double swsl2 = scaleweight2[sl];
                double kbasesl = exp(-v/rsl2);
                double ksl = Ks(kbasesl,swsl2);

                double aislTXalsi; VECDOT(aislTXalsi,aisl,alsl);

                Ett = Ett+ksl*aislTXalsi;
            }
        }
    }
    } else {
#pragma omp parallel for schedule(static) shared(rhot,L,R,dim,CSP,CSPL,scales2,scaleweight2) private(i,l,sl) reduction(+:Ett)
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

                double aislTXalsi; _3VECDOT(aislTXalsi,aisl,alsl);

                Ett = Ett+ksl*aislTXalsi;
            }
        }
    }
    }

    Et[0] = Ett;
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *rhot;
    int L;
    int R;
    int dim;
    int CSP;
    double *scales2;
    double *scaleweight2;
    double energyweight;
    double *Et;

    /* check for proper number of arguments */
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("gradC:nrhs","6 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("gradC:nlhs","One output required.");
    }
    
    /* get the values of the non-array inputs  */
    L = (double)(mxGetScalar(prhs[1]));
    R = (double)(mxGetScalar(prhs[2]));
    dim = (double)(mxGetScalar(prhs[3]));

    CSP = dim*(1+R);

    /* get the array inputs */
    rhot = mxGetPr(prhs[0]);
    if(CSP*L != mxGetM(prhs[0])) {
        mexErrMsgIdAndTxt("gradC:nlhs","rhot dimension mismatch (got %d, expected %d).",mxGetM(prhs[0]),CSP*L);
    }
    scales2 = mxGetPr(prhs[4]);
    if(R != mxGetN(prhs[4])) {
        mexErrMsgIdAndTxt("gradC:nlhs","scales2 dimension mismatch (got %d, expected %d).",mxGetN(prhs[4]),R);
    }
    scaleweight2 = mxGetPr(prhs[5]);
    if(R != mxGetN(prhs[5])) {
        mexErrMsgIdAndTxt("gradC:nlhs","scaleweight2 dimension mismatch (got %d, expected %d).",mxGetN(prhs[5]),R);
    }

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    Et = mxGetPr(plhs[0]);

    /* call the computational routine */
    energyScaleC(Et,rhot,L,R,dim,CSP,scales2,scaleweight2);
}
