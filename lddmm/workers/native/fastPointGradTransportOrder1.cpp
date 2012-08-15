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

#include "fastPointGradTransportOrder1.hpp"
#include "dkernelDerivativeMatrix.hpp"

using namespace sla;

//void addEgradt(scalar *dy, scalar *yt, scalar *Gt, scalar *rhoj0, int L, int R, int CSP, int CSD, scalar *scales2, scalar *scaleweight2, scalar energyweight)
//{
//    int CSPL = CSP*L;
//    const int dim = 3;
//
//    int i, l, sl, si;
//
//}

/* The computational routine */
void fastPointGradTransportOrder1(scalar *dy, scalar *yt, scalar *Gt, scalar *rhoj0, int L, int R, scalar *scales2, scalar *scaleweight2, scalar energyweight)
{
    /* dy already initialized to zeros */

#pragma omp parallel for schedule(static) shared(dy,yt,Gt,rhoj0,L,R,scales2,scaleweight2,energyweight)
    for (int i=0; i<L; i++) { /* grid point */

        Vector3<scalar> phii;
        Matrix3<scalar> DphiiT;
        Vector3<scalar> mui;
        Matrix3<scalar> mujiT;

        for (int n=0; n<L; n++) { /* particle */
            Vector3<scalar> dphin(&yt[INDPHI(n)]);
            Vector3<scalar> dmun(&yt[INDMU(n)]);

            // phi parts
            phii = phii + aphiphiT(n,i,dphin,L,Gt,rhoj0,scales2[0],scaleweight2[0]) + amuphiT(n,i,dmun,L,Gt,rhoj0,scales2[0],scaleweight2[0]);

            for (int l=0; l<dim; l++) {
                Vector3<scalar> ddphinl(&yt[INDDPHI(n)+l*dim]);

                phii = phii + aDphilphiT(n,i,l,ddphinl,L,Gt,rhoj0,scales2[0],scaleweight2[0]);
            }

            for (int j=0; j<dim; j++) {
                Vector3<scalar> dmujnj(&yt[INDMUJ(n)+j*dim]);

                phii = phii + amujphiT(n,i,j,dmujnj,L,Gt,rhoj0,scales2[0],scaleweight2[0]);
            }

            for (int k=0; k<dim; k++) {
                DphiiT.addRow(k,
                        aphiDphilT(n,i,k,dphin,L,Gt,rhoj0,scales2[0],scaleweight2[0])
                        +amuDphilT(n,i,k,dmun,L,Gt,rhoj0,scales2[0],scaleweight2[0])
                        );

                for (int l=0; l<dim; l++) {
                    Vector3<scalar> ddphinl(&yt[INDDPHI(n)+l*dim]);

                    DphiiT.addRow(k,
                            aDphilDphikT(n,i,l,k,ddphinl,L,Gt,rhoj0,scales2[0],scaleweight2[0])
                            );
                }

                for (int j=0; j<dim; j++) {
                    Vector3<scalar> dmujnj(&yt[INDMUJ(n)+j*dim]);

                    DphiiT.addRow(k,
                            amujDphilT(n,i,j,k,dmujnj,L,Gt,rhoj0,scales2[0],scaleweight2[0])
                            );
                }

            }


            // mu parts
            mui = mui + aphimuT(n,i,dphin,L,Gt,rhoj0,scales2[0],scaleweight2[0]) + amumuT(n,i,dmun,L,Gt,rhoj0,scales2[0],scaleweight2[0]);

            for (int l=0; l<dim; l++) {
                Vector3<scalar> ddphinl(&yt[INDDPHI(n)+l*dim]);

                mui = mui + aDphilmuT(n,i,l,ddphinl,L,Gt,rhoj0,scales2[0],scaleweight2[0]);
            }

            for (int j=0; j<dim; j++) {
                Vector3<scalar> dmujnj(&yt[INDMUJ(n)+j*dim]);

                mui = mui + amujmuT(n,i,j,dmujnj,L,Gt,rhoj0,scales2[0],scaleweight2[0]);
            }

            for (int jm=0; jm<dim; jm++) {
                mujiT.addRow(jm,
                        aphimujT(n,i,jm,dphin,L,Gt,rhoj0,scales2[0],scaleweight2[0])
                        +amumujT(n,i,jm,dmun,L,Gt,rhoj0,scales2[0],scaleweight2[0])
                        );

                for (int l=0; l<dim; l++) {
                    Vector3<scalar> ddphinl(&yt[INDDPHI(n)+l*dim]);

                    mujiT.addRow(jm,
                            aDphilmujT(n,i,l,jm,ddphinl,L,Gt,rhoj0,scales2[0],scaleweight2[0])
                            );
                }

                for (int j=0; j<dim; j++) {
                    Vector3<scalar> dmujnj(&yt[INDMUJ(n)+j*dim]);

                    mujiT.addRow(jm,
                            amujmujT(n,i,j,jm,dmujnj,L,Gt,rhoj0,scales2[0],scaleweight2[0])
                            );
                }

            }

        }

        phii.store(&dy[INDPHI(i)]);
        DphiiT.store(&dy[INDDPHI(i)]);
        mui.store(&dy[INDMU(i)]);
        mujiT.store(&dy[INDMUJ(i)]);
    }

//    addEgradt(dy,yt,rhot,L,R,CSD,scales2,scaleweight2,energyweight);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *yt;
    double *Gt;
    double *rhoj0;
    int L;
    int R;
    int ldim;
    double *scales2;
    double *scaleweight2;
    double energyweight;
    double *dy;

    /* check for proper number of arguments */
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("gradC:nrhs","9 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("gradC:nlhs","One output required.");
    }
    
    /* get the values of the non-array inputs  */
    L = (double)(mxGetScalar(prhs[3]));
    R = (double)(mxGetScalar(prhs[4]));
    ldim = (double)(mxGetScalar(prhs[5]));
    mexAssert(ldim == dim);
    energyweight = (double)(mxGetScalar(prhs[8]));

    /* get the array inputs */
    yt = mxGetPr(prhs[0]);
    int n= mxGetM(prhs[0]);
    if(n != gradCSP*L) {
        mexErrMsgIdAndTxt("gradC:nlhs","yt dimension mismatch (got %d, expected %d).",mxGetM(prhs[0]),gradCSP*L);
    }
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
    dy = mxGetPr(plhs[0]);

    /* debug */
    /*mexPrintf("%f %f %f %d %d %d %d %f %f %d %d %d\n",t,yt[0],rhot[0],L,R,CSP,CSD,scales2[0],scaleweight2[0],GL,GS,GI);*/

    /* call the computational routine */
    fastPointGradTransportOrder1(dy,yt,Gt,rhoj0,L,R,scales2,scaleweight2,energyweight);
}
