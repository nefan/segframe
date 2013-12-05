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

#pragma omp parallel for schedule(static) shared(K__ij,D1Ks__ijb,D2Ks__ijbg,D3Ks__ijbgd)
    for (int j=0; j<Lp; j++) {
        Vector3<scalar> xj(&p_a_i.x[p_a_i.rows*j]);

        for (int i=0; i<Lq; i++) {
            Vector3<scalar> xi(&q_a_i.x[q_a_i.rows*i]);
            Vector3<scalar> ximxj = xi-xj;

            const scalar r = sqrt(scales2.x[sl]);

            if (sqrt(dot(ximxj,ximxj)) > 4*r) // cutoff
                continue;

            const scalar ks = gs::gamma(ximxj,scales2.x[sl],scaleweight2.x[sl]);
            K__ij.x[i+K__ij.rows*j] = ks;

            if (nargout > 1) {

                for (int b=0; b<dim; b++) {
                    Vector3<int> da; da.set(1,b);
                    D1Ks__ijb.x[i+D1Ks__ijb.dimsI[0]*j+D1Ks__ijb.dimsI[1]*b] =
                        gs::DaKs(da,ximxj,r,ks);

                    if (nargout > 2) {

                        for (int g=0; g<dim; g++) {
                            Vector3<int> db = da; db.set(db[g]+1,g);
                            D2Ks__ijbg.x[i+D2Ks__ijbg.dimsI[0]*j+D2Ks__ijbg.dimsI[1]*b+D2Ks__ijbg.dimsI[2]*g] =
                                gs::DaKs(db,ximxj,r,ks);

                            if (nargout > 3) {

                                for (int d=0; d<dim; d++) {
                                    Vector3<int> dc = db; dc.set(dc[d]+1,d);
                                    D3Ks__ijbgd.x[i+D3Ks__ijbgd.dimsI[0]*j+D3Ks__ijbgd.dimsI[1]*b+D3Ks__ijbgd.dimsI[2]*g+D3Ks__ijbgd.dimsI[3]*d] =
                                        gs::DaKs(dc,ximxj,r,ks);

                                    if (nargout > 4) {

                                        for (int e=0; e<dim; e++) {
                                            Vector3<int> dd = dc; dd.set(dd[e]+1,e);
                                            D4Ks__ijbgde.x[i+D4Ks__ijbgde.dimsI[0]*j+D4Ks__ijbgde.dimsI[1]*b+D4Ks__ijbgde.dimsI[2]*g+D4Ks__ijbgde.dimsI[3]*d+D4Ks__ijbgde.dimsI[4]*e] =
                                                gs::DaKs(dd,ximxj,r,ks);

                                            if (nargout > 5) {

                                                for (int phi=0; phi<dim; phi++) {
                                                    Vector3<int> de = dd; de.set(de[phi]+1,phi);
                                                    D5Ks__ijbgdep.x[i+D5Ks__ijbgdep.dimsI[0]*j+D5Ks__ijbgdep.dimsI[1]*b+D5Ks__ijbgdep.dimsI[2]*g+D5Ks__ijbgdep.dimsI[3]*d+D5Ks__ijbgdep.dimsI[4]*e+D5Ks__ijbgdep.dimsI[5]*phi] =
                                                        gs::DaKs(de,ximxj,r,ks);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

cljmex_end()
