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
void scalesJac(double *jac, double *grid, double *rho, int scaleI, int L, int R, int dim, int CSP, double *scales2, double *scaleweight2, int Ngrid)
{
    /* jac already initialized to zeros */
    int g, i, j, l;

    if (dim == 2) {
        /* not implemented */
        mexAssert(false);
    } else { /* dim == 3 */
#pragma omp parallel for schedule(static) shared(jac,grid,rho,scaleI,L,R,dim,CSP,scales2,scaleweight2) private(g,i,j,l)
    for (g=0; g<Ngrid; g++) { /* grid point */
        double x[3] = {grid[dim*g+0],grid[dim*g+1],grid[dim*g+2]};

        for (i=0; i<dim; i++) {
            for (j=0; j<dim; j++) {

                double v = 0;
                for (l=0; l<L; l++) { /* particle */
                    double xl[3] = {rho[INDRHOX(l,0)],rho[INDRHOX(l,1)],rho[INDRHOX(l,2)]};
                    double ali = rho[INDRHOP(l,i,scaleI)];

                    double xmxl[3]; _3VECMINUS(xmxl,x,xl);
                    double w = _3VECDOT2(xmxl,xmxl);

                    double r2 = scales2[scaleI];
                    double sw2 = scaleweight2[scaleI];
                    double kbase = exp(-w/r2);
                    double d1ks = D1Ks(kbase,r2,sw2);

                    v += ali*d1ks*2*xmxl[j];
                }

                jac[g] += v*v;
            }
        }
    }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *grid;
    double *rho;
    int scaleI;
    int L;
    int R;
    int dim;
    double *scales2;
    double *scaleweight2;
    double *jac;

    /* check for proper number of arguments */
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("scalesJac:nrhs","8 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("scalesJac:nlhs","One output required.");
    }
    
    /* get the values of the non-array inputs  */
    scaleI = (double)(mxGetScalar(prhs[2]))-1;
    L = (double)(mxGetScalar(prhs[3]));
    R = (double)(mxGetScalar(prhs[4]));
    dim = (double)(mxGetScalar(prhs[5]));

    int CSP = dim*(1+R);

    /* get the array inputs */
    grid = mxGetPr(prhs[0]);
    int n= mxGetM(prhs[0]);
    int Ngrid = n/3;
    rho = mxGetPr(prhs[1]);
    if(CSP*L != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("gradC:nlhs","rho dimension mismatch (got %d, expected %d).",mxGetM(prhs[1]),CSP*L);
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
    plhs[0] = mxCreateDoubleMatrix(Ngrid,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    jac = mxGetPr(plhs[0]);

    /* debug */
    /*mexPrintf("%f %f %f %d %d %d %d %f %f %d %d %d\n",t,yt[0],rhot[0],L,R,CSP,CSD,scales2[0],scaleweight2[0],GL,GS,GI);*/

    /* call the computational routine */
    scalesJac(jac,grid,rho,scaleI,L,R,dim,CSP,scales2,scaleweight2,Ngrid);
}
