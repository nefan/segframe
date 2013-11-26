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

//#include "gaussianTKsC.hpp"
#include "gaussianTKsC.cljmex.hpp"
#include "common.h"
#include "gaussianKernels.hpp"

/* 
 * arrays of Gaussian kernels with derivatives
 */
cljmex_start()
    assert(dim == 3);
    assert(R == 1);
    const int sl = 0;
    assert(order == 0);

    // hardcoded - move to cljmex
    mwSize D1Ksdim[3] = {L, L, dim};
    mwSize D1KsdimI[2] = {D1Ksdim[0], D1Ksdim[0]*D1Ksdim[1]};
    pargout[1] = (mxArray *)mxCreateNumericArray(3,D1Ksdim,mxDOUBLE_CLASS,mxREAL);
    double *D1Ks__ijb = mxGetPr (pargout[1]);

    mwSize D2Ksdim[4] = {L, L, dim, dim};
    mwSize D2KsdimI[3] = {D2Ksdim[0], D2Ksdim[0]*D2Ksdim[1], D2Ksdim[0]*D2Ksdim[1]*D2Ksdim[2]};
    pargout[2] = (mxArray *)mxCreateNumericArray(4,D2Ksdim,mxDOUBLE_CLASS,mxREAL);
    double *D2Ks__ijbg = mxGetPr (pargout[2]);

    mwSize D3Ksdim[5] = {L, L, dim, dim, dim};
    mwSize D3KsdimI[4] = {D3Ksdim[0], D3Ksdim[0]*D3Ksdim[1], D3Ksdim[0]*D3Ksdim[1]*D3Ksdim[2], D3Ksdim[0]*D3Ksdim[1]*D3Ksdim[2]*D3Ksdim[3]};
    pargout[3] = (mxArray *)mxCreateNumericArray(5,D3Ksdim,mxDOUBLE_CLASS,mxREAL);
    double *D3Ks__ijbgd = mxGetPr (pargout[3]);

    // sum the entries in the sparse array R
    for (int j=0; j<L; j++)
        for (int i=0; i<L; i++) {
            Vector3<scalar> xi(&q0_a_i.x[q0_a_i.rows*i]);
            Vector3<scalar> xj(&q0_a_i.x[q0_a_i.rows*j]);

            const scalar ks = gs::gamma(xi-xj,scales2.x[sl],scaleweight2.x[sl]);
            K__ij.x[i+K__ij.rows*j] = ks;

            if (nargout > 1) {

                for (int b=0; b<dim; b++) {
                    Vector3<int> da; da.set(1,b);
                    D1Ks__ijb[i+D1KsdimI[0]*j+D1KsdimI[1]*b] =
                        gs::DaKs(da,xi-xj,sqrt(scales2.x[sl]),ks);

                    if (nargout > 2) {

                        for (int g=0; g<dim; g++) {
                            Vector3<int> db = da; db.set(db[g]+1,g);
                            D2Ks__ijbg[i+D2KsdimI[0]*j+D2KsdimI[1]*b+D2KsdimI[2]*g] =
                                gs::DaKs(db,xi-xj,sqrt(scales2.x[sl]),ks);

                            if (nargout > 3) {

                                for (int d=0; d<dim; d++) {
                                    Vector3<int> dc = db; dc.set(dc[d]+1,d);
                                    D3Ks__ijbgd[i+D3KsdimI[0]*j+D3KsdimI[1]*b+D3KsdimI[2]*g+D3KsdimI[3]*d] =
                                        gs::DaKs(dc,xi-xj,sqrt(scales2.x[sl]),ks);
                                }
                            }

                        }
                    }
                }
            }
        }

//       plhs[0] = (mxArray *) mxCreateNumericArray 
//         (mxGetNumberOfDimensions (prhs[0]),
//          mxGetDimensions (prhs[0]), mxGetClassID (prhs[0]),
//          mxIsComplex (prhs[0]));
//       vri = mxGetPr (prhs[0]);

cljmex_end()
