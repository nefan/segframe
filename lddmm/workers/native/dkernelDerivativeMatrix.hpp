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

using namespace sla;
using namespace gs;
using namespace kd;

// indexing
inline Vector3<scalar> x(const int i, const scalar* Gt) {
    Vector3<scalar> xi(Gt[INDGX(i,0)],Gt[INDGX(i,1)],Gt[INDGX(i,2)]);
    
    return xi;
}
inline Vector3<scalar> mu(const int i, const scalar* Gt) {
    Vector3<scalar> mui(Gt[INDGMU(i,0)],Gt[INDGMU(i,1)],Gt[INDGMU(i,2)]);
    
    return mui;
}
inline Matrix3<scalar> DphiT(const int i, const scalar* Gt) {
    Matrix3<scalar> DphiiT(
        Vector3<scalar>(Gt[INDGDPHI(i,0,0)],Gt[INDGDPHI(i,1,0)],Gt[INDGDPHI(i,2,0)]),
        Vector3<scalar>(Gt[INDGDPHI(i,0,1)],Gt[INDGDPHI(i,1,1)],Gt[INDGDPHI(i,2,1)]),
        Vector3<scalar>(Gt[INDGDPHI(i,0,2)],Gt[INDGDPHI(i,1,2)],Gt[INDGDPHI(i,2,2)]));
    
    return DphiiT;
}
inline Matrix3<scalar> Dphi(const int i, const scalar* Gt) {
    
    return Dphi(i,Gt).T();
}
inline Vector3<scalar> Dphil(const int i, const int l, const scalar* Gt) {
    Vector3<scalar> dphiil(Gt[INDGDPHI(i,0,l)],Gt[INDGDPHI(i,1,l)],Gt[INDGDPHI(i,2,l)]);
    
    return dphiil;
}
inline Vector3<scalar> muj(const int i, const int j, const scalar* Gt, const scalar* rhoj0) {
    Vector3<scalar> rhoj0j(rhoj0[INDRHOJ(i,0,j)],rhoj0[INDRHOJ(i,1,j)],rhoj0[INDRHOJ(i,2,j)]);

    Vector3<scalar> muij = DphiT(i,Gt).inv()*rhoj0j;
    
    return muij;
}
inline Vector3<scalar> dtDphi(const int n, const Vector3<scalar>& v, const scalar* Gt, const scalar* rhoj0, const int L, const scalar r2, const scalar sw2) {
    Vector3<scalar> xn = x(n,Gt);

    Vector3<scalar> a;

    for (int l=0; l<dim; l++) {
        Vector3<scalar> Dphinl = Dphil(n,l,Gt);

        for (int i=0; i<L; i++) {
            Vector3<scalar> xi = x(i,Gt);

            Vector3<scalar> ximxn = xi-xn;

            a = a + v[l]*dot<scalar>(Dphinl,N2Ks(ximxn,r2,sw2))*mu(i,Gt);
            for (int j=0; j<dim; j++) {
                Vector3<scalar> Dphiij = Dphil(i,j,Gt);

                a = a + v[l]*dot<scalar>(Dphinl,D1N2Ks(ximxn,Dphiij,r2,sw2))*muj(i,j,Gt,rhoj0);
            }
        }
    }
    
    return a;
}
inline Vector3<scalar> dtDphiT(const int n, const Vector3<scalar>& v, const scalar* Gt, const scalar* rhoj0, const int L, const scalar r2, const scalar sw2) {
    Vector3<scalar> xn = x(n,Gt);

    Vector3<scalar> a;

    for (int l=0; l<dim; l++) {
        Vector3<scalar> Dphinl = Dphil(n,l,Gt);
        scalar al = 0;

        for (int i=0; i<L; i++) {
            Vector3<scalar> xi = x(i,Gt);

            Vector3<scalar> ximxn = xi-xn;

            al += dot<scalar>(Dphinl,N2Ks(ximxn,r2,sw2))*dot<scalar>(mu(i,Gt),v);
            for (int j=0; j<dim; j++) {
                Vector3<scalar> Dphiij = Dphil(i,j,Gt);

                al += dot<scalar>(Dphinl,D1N2Ks(ximxn,Dphiij,r2,sw2))*dot<scalar>(muj(i,j,Gt,rhoj0),v);
            }
        }

        a.set(al,l);
    }
    
    return a;
}

inline Vector3<scalar> aphiphiT(const int n, const int i, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> mui = mu(i,Gt);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = (2.0*d1gamma(ximxn,r2,sw2)*dot<scalar>(mui,v))*ximxn;

    for (int jm = 0; jm < dim; jm++) {
        Vector3<scalar> Dphiijm = Dphil(i,jm,Gt);
        a = a + dot<scalar>(muj(i,jm,Gt,rhoj0),v)*(
                2.0*d1gamma(ximxn,r2,sw2)*Dphiijm
                +4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphiijm,ximxn)*ximxn);

    }

    if (i == n) {
        for (int im = 0; im<L; im++) {
            Vector3<scalar> xim = x(im,Gt);
            Vector3<scalar> muim = mu(im,Gt);

            Vector3<scalar> ximmxn = xim-xn;

            a = a - (2.0*d1gamma(ximmxn,r2,sw2)*dot<scalar>(muim,v))*ximmxn;

            for (int jm = 0; jm < dim; jm++) {
                Vector3<scalar> Dphiimjm = Dphil(im,jm,Gt);
                a = a - dot<scalar>(muj(im,jm,Gt,rhoj0),v)*(
                        2.0*d1gamma(ximmxn,r2,sw2)*Dphiimjm
                        +4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(Dphiimjm,ximmxn)*ximmxn);

            }
        }
    }

    return a;
}

inline Vector3<scalar> aphiDphilT(const int n, const int i, const int l, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> mujil = muj(i,l,Gt,rhoj0);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = dot<scalar>(mujil,v)*N1Ks(ximxn,r2,sw2);

    return a;
}

inline Vector3<scalar> aphimuT(const int n, const int i, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = gamma(ximxn,r2,sw2)*v;

    return a;
}

inline Vector3<scalar> aphimujT(const int n, const int i, const int j, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> Dphiij = Dphil(i,j,Gt);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = dot<scalar>(Dphiij,N1Ks(ximxn,r2,sw2))*v;

    return a;
}

inline Vector3<scalar> aDphilphiT(const int n, const int i, const int l, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> mui = mu(i,Gt);
    Vector3<scalar> Dphinl = Dphil(n,l,Gt);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = -dot<scalar>(mui,v)*((4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,Dphinl))*ximxn+(2.0*d1gamma(ximxn,r2,sw2))*Dphinl);

    for (int jm = 0; jm < dim; jm++) {
        Vector3<scalar> Dphiijm = Dphil(i,jm,Gt);
        a = a - dot<scalar>(muj(i,jm,Gt,rhoj0),v)*(
                +4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphiijm,Dphinl)*ximxn
                +8.0*d3gamma(ximxn,r2,sw2)*dot<scalar>(Dphiijm,ximxn)*dot<scalar>(Dphinl,ximxn)*ximxn
                +4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphiijm,ximxn)*Dphinl
                +4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphinl,ximxn)*Dphiijm);

    }

    if (i == n) {
        for (int im = 0; im<L; im++) {
            Vector3<scalar> xim = x(im,Gt);
            Vector3<scalar> muim = mu(im,Gt);

            Vector3<scalar> ximmxn = xim-xn;

            a = a + dot<scalar>(muim,v)*((4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(ximmxn,Dphinl))*ximmxn+(2.0*d1gamma(ximmxn,r2,sw2))*Dphinl);

            for (int jm = 0; jm < dim; jm++) {
                Vector3<scalar> Dphiimjm = Dphil(im,jm,Gt);

                a = a + dot<scalar>(muj(im,jm,Gt,rhoj0),v)*(
                        +4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(Dphiimjm,Dphinl)*ximmxn
                        +8.0*d3gamma(ximmxn,r2,sw2)*dot<scalar>(Dphiimjm,ximmxn)*dot<scalar>(Dphinl,ximmxn)*ximmxn
                        +4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(Dphiimjm,ximmxn)*Dphinl
                        +4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(Dphinl,ximmxn)*Dphiimjm);
            }
        }
    }

    return a;
}

inline Vector3<scalar> aDphilDphikT(const int n, const int i, const int l, const int k, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> mujik = muj(i,k,Gt,rhoj0);
    Vector3<scalar> Dphinl = Dphil(n,l,Gt);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = dot<scalar>(mujik,v)*D1N2KsT(ximxn,Dphinl,r2,sw2);

    if (i == n && l == k) {
        for (int im = 0; im < L; im++) {
            Vector3<scalar> xim = x(im,Gt);
            Vector3<scalar> ximmxn = xim-xn;

            a = a + dot<scalar>(mu(im,Gt),v)*N2Ks(ximmxn,r2,sw2);

            for (int jm = 0; jm < dim; jm++) {
                Vector3<scalar> Dphiimjm = Dphil(im,jm,Gt);

                a = a + dot<scalar>(muj(im,jm,Gt,rhoj0),v)*D1N2Ks(ximmxn,Dphiimjm,r2,sw2);
            }
        }
    }

    return a;
}

inline Vector3<scalar> aDphilmuT(const int n, const int i, const int l, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> Dphinl = Dphil(n,l,Gt);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = dot<scalar>(Dphinl,N2Ks(ximxn,r2,sw2))*v;

    return a;
}

inline Vector3<scalar> aDphilmujT(const int n, const int i, const int l, const int j, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> Dphinl = Dphil(n,l,Gt);
    Vector3<scalar> Dphiij = Dphil(i,j,Gt);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = dot<scalar>(Dphinl,D1N2Ks(ximxn,Dphiij,r2,sw2))*v;

    return a;
}

inline Vector3<scalar> amuphiT(const int n, const int i, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> mui = mu(i,Gt);
    Vector3<scalar> mun = mu(n,Gt);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = dot<scalar>(mun,mui)*(
            2.0*d1gamma(ximxn,r2,sw2)*v
            +4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,v)*ximxn);

    for (int jm = 0; jm < dim; jm++) {
        Vector3<scalar> Dphinjm = Dphil(n,jm,Gt);
        Vector3<scalar> mujnjm = muj(n,jm,Gt,rhoj0);
        Vector3<scalar> mujijm = muj(i,jm,Gt,rhoj0);

        a = a - (dot<scalar>(mujnjm,mui)-dot<scalar>(mun,mujijm))*(
                4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphinjm,v)*ximxn
                +8.0*d3gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,v)*dot<scalar>(Dphinjm,ximxn)*ximxn
                +4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphinjm,ximxn)*v
                +4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,v)*Dphinjm);

        for (int j = 0; j < dim; j++) {
            Vector3<scalar> Dphiij = Dphil(i,j,Gt);
            Vector3<scalar> mujij = muj(i,j,Gt,rhoj0);

            a = a - dot<scalar>(mujnjm,mujij)*(
                    +8.0*d3gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,Dphinjm)*dot<scalar>(Dphiij,v)*ximxn
                    +8.0*d3gamma(ximxn,r2,sw2)*dot<scalar>(Dphiij,Dphinjm)*dot<scalar>(ximxn,v)*ximxn
                    +4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphiij,v)*Dphinjm
                    +4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphiij,Dphinjm)*v
                    +8.0*d3gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,Dphiij)*dot<scalar>(Dphinjm,v)*ximxn
                    +4.0*d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphinjm,v)*Dphiij
                    +16.0*d4gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,Dphiij)*dot<scalar>(ximxn,Dphinjm)*dot<scalar>(ximxn,v)*ximxn
                    +8.0*d3gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,Dphinjm)*dot<scalar>(ximxn,v)*Dphiij
                    +8.0*d3gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,Dphiij)*dot<scalar>(ximxn,Dphinjm)*v
                    +8.0*d3gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,Dphiij)*dot<scalar>(ximxn,v)*Dphinjm
                    );
        }
    }

    if (i == n) {
        for (int im = 0; im < L; im++) {
            Vector3<scalar> xim = x(im,Gt);
            Vector3<scalar> muim = mu(im,Gt);

            Vector3<scalar> ximmxn = xim-xn;

            a = a - dot<scalar>(mun,muim)*(
                    2.0*d1gamma(ximmxn,r2,sw2)*v
                    +4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(ximmxn,v)*ximmxn);

            for (int jm = 0; jm < dim; jm++) {
                Vector3<scalar> Dphinjm = Dphil(n,jm,Gt);
                Vector3<scalar> mujnjm = muj(n,jm,Gt,rhoj0);
                Vector3<scalar> mujimjm = muj(im,jm,Gt,rhoj0);

                a = a + (dot<scalar>(mujnjm,muim)-dot<scalar>(mun,mujimjm))*(
                        4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(Dphinjm,v)*ximmxn
                        +8.0*d3gamma(ximmxn,r2,sw2)*dot<scalar>(ximmxn,v)*dot<scalar>(Dphinjm,ximmxn)*ximmxn
                        +4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(Dphinjm,ximmxn)*v
                        +4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(ximmxn,v)*Dphinjm);

                for (int j = 0; j < dim; j++) {
                    Vector3<scalar> Dphiimj = Dphil(im,j,Gt);
                    Vector3<scalar> mujimj = muj(im,j,Gt,rhoj0);

                    a = a + dot<scalar>(mujnjm,mujimj)*(
                            +8.0*d3gamma(ximmxn,r2,sw2)*dot<scalar>(ximmxn,Dphinjm)*dot<scalar>(Dphiimj,v)*ximmxn
                            +8.0*d3gamma(ximmxn,r2,sw2)*dot<scalar>(Dphiimj,Dphinjm)*dot<scalar>(ximmxn,v)*ximmxn
                            +4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(Dphiimj,v)*Dphinjm
                            +4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(Dphiimj,Dphinjm)*v
                            +8.0*d3gamma(ximmxn,r2,sw2)*dot<scalar>(ximmxn,Dphiimj)*dot<scalar>(Dphinjm,v)*ximmxn
                            +4.0*d2gamma(ximmxn,r2,sw2)*dot<scalar>(Dphinjm,v)*Dphiimj
                            +16.0*d4gamma(ximmxn,r2,sw2)*dot<scalar>(ximmxn,Dphiimj)*dot<scalar>(ximmxn,Dphinjm)*dot<scalar>(ximmxn,v)*ximmxn
                            +8.0*d3gamma(ximmxn,r2,sw2)*dot<scalar>(ximmxn,Dphinjm)*dot<scalar>(ximmxn,v)*Dphiimj
                            +8.0*d3gamma(ximmxn,r2,sw2)*dot<scalar>(ximmxn,Dphiimj)*dot<scalar>(ximmxn,Dphinjm)*v
                            +8.0*d3gamma(ximmxn,r2,sw2)*dot<scalar>(ximmxn,Dphiimj)*dot<scalar>(ximmxn,v)*Dphinjm
                            );
                }
            }
        }
    }

    return a;
}

inline Vector3<scalar> amuDphilT(const int n, const int i, const int l, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> mujil = muj(i,l,Gt,rhoj0);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a;

    for (int jm = 0; jm < dim; jm++) {
        Vector3<scalar> Dphinjm = Dphil(n,jm,Gt);

        a = a - 4.0*dot<scalar>(muj(n,jm,Gt,rhoj0),mujil)*(
                +d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphinjm,ximxn)*v
                +d2gamma(ximxn,r2,sw2)*dot<scalar>(ximxn,v)*Dphinjm
                +d2gamma(ximxn,r2,sw2)*dot<scalar>(Dphinjm,v)*ximxn
                +2.0*d3gamma(ximxn,r2,sw2)*dot<scalar>(Dphinjm,ximxn)*dot<scalar>(ximxn,v)*ximxn);
    }

    if (i == n) {
        for (int im = 0; im < L; im++) {
            Vector3<scalar> xim = x(im,Gt);
            Vector3<scalar> ximmxn = xim-xn;

            a = a - (dot<scalar>(muj(n,l,Gt,rhoj0),mu(im,Gt))-dot<scalar>(mu(n,Gt),muj(im,l,Gt,rhoj0)))
                *D2N2KsT(ximmxn,v,r2,sw2);

            for (int jm = 0; jm < dim; jm++) {
                Vector3<scalar> Dphiimjm = Dphil(im,jm,Gt);

                a = a - dot<scalar>(muj(n,l,Gt,rhoj0),muj(im,jm,Gt,rhoj0))
                *D2D1N2KsT(ximmxn,Dphiimjm,v,r2,sw2);
            }
        }
    }

    return a;
}

inline Vector3<scalar> amumuT(const int n, const int i, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> mun = mu(n,Gt);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = -dot<scalar>(v,N2Ks(ximxn,r2,sw2))*mun;

    for (int jm = 0; jm < dim; jm++) {
        Vector3<scalar> Dphinjm = Dphil(n,jm,Gt);
        a = a - dot<scalar>(D2N2Ks(ximxn,Dphinjm,r2,sw2),v)*muj(n,jm,Gt,rhoj0);

    }

    if (i == n) {
        for (int im = 0; im<L; im++) {
            Vector3<scalar> xim = x(im,Gt);
            Vector3<scalar> ximmxn = xim-xn;

            a = a - dot<scalar>(v,N2Ks(ximmxn,r2,sw2))*mu(im,Gt);

            for (int jm = 0; jm < dim; jm++) {
                Vector3<scalar> Dphinjm = Dphil(n,jm,Gt);
                a = a + dot<scalar>(D2N2Ks(ximmxn,Dphinjm,r2,sw2),v)*muj(im,jm,Gt,rhoj0);

            }
        }
    }

    return a;
}

inline Vector3<scalar> amumujT(const int n, const int i, const int j, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> xi = x(i,Gt);
    Vector3<scalar> xn = x(n,Gt);
    Vector3<scalar> mun = mu(n,Gt);
    Vector3<scalar> Dphiij = Dphil(i,j,Gt);
    Vector3<scalar> Dphinj = Dphil(n,j,Gt);

    Vector3<scalar> ximxn = xi-xn;

    Vector3<scalar> a = dot<scalar>(v,D2N2Ks(ximxn,Dphinj,r2,sw2))*mun;

    for (int jm = 0; jm < dim; jm++) {
        Vector3<scalar> Dphinjm = Dphil(n,jm,Gt);
        a = a - dot<scalar>(D2D1N2Ks(ximxn,Dphiij,Dphinjm,r2,sw2),v)*muj(n,jm,Gt,rhoj0);

    }

    if (i == n) {
        for (int im = 0; im<L; im++) {
            Vector3<scalar> xim = x(im,Gt);
            Vector3<scalar> ximmxn = xim-xn;
            Vector3<scalar> muim = mu(im,Gt);

            a = a - dot<scalar>(v,D2N2Ks(ximmxn,Dphinj,r2,sw2))*muim;

            for (int jm = 0; jm < dim; jm++) {
                Vector3<scalar> Dphiimjm = Dphil(im,jm,Gt);
                a = a - dot<scalar>(D2D1N2Ks(ximmxn,Dphiimjm,Dphinj,r2,sw2),v)*muj(im,jm,Gt,rhoj0);

            }

        }
    }

    return a;
}

// obvious optimization: calls matrix
inline Vector3<scalar> amujphiT(const int n, const int i, const int j, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> mujnj = muj(n,j,Gt,rhoj0);
    Matrix3<scalar> DphiiTinv = DphiT(i,Gt).inv();

    Vector3<scalar> a;

    for (int jm = 0; jm < dim; jm++) {
        a = a - dot<scalar>(DphiiTinv.column(jm),v)*aDphilphiT(n,i,jm,mujnj,L,Gt,rhoj0,r2,sw2);

    }

    return a;
}

// obvious optimization: calls matrix
inline Vector3<scalar> amujDphilT(const int n, const int i, const int j, const int l, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> mujnj = muj(n,j,Gt,rhoj0);
    Matrix3<scalar> DphiiTinv = DphiT(i,Gt).inv();

    Vector3<scalar> a;

    for (int jm = 0; jm < dim; jm++) {
        a = a - dot<scalar>(DphiiTinv.column(jm),v)*aDphilDphikT(n,i,jm,l,mujnj,L,Gt,rhoj0,r2,sw2);

    }

    if (i == n) {
        // DphiiTinv == DphinTinv
        a = a + dot<scalar>(DphiiTinv.column(l),v)*(DphiiTinv*dtDphiT(n,mujnj,Gt,rhoj0,L,r2,sw2));
    }

    return a;
}

// obvious optimization: calls matrix
inline Vector3<scalar> amujmuT(const int n, const int i, const int j, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> mujnj = muj(n,j,Gt,rhoj0);
    Matrix3<scalar> DphiiTinv = DphiT(i,Gt).inv();

    Vector3<scalar> a;

    for (int jm = 0; jm < dim; jm++) {
        a = a - dot<scalar>(DphiiTinv.column(jm),v)*aDphilmuT(n,i,jm,mujnj,L,Gt,rhoj0,r2,sw2);

    }

    return a;
}

// obvious optimization: calls matrix
inline Vector3<scalar> amujmujT(const int n, const int i, const int j, const int jm, const Vector3<scalar> &v, const int L, const scalar* Gt, const scalar* rhoj0, const scalar r2, const scalar sw2) {
    Vector3<scalar> mujnj = muj(n,j,Gt,rhoj0);
    Matrix3<scalar> DphiiTinv = DphiT(i,Gt).inv();

    Vector3<scalar> a;

    for (int jmm = 0; jmm < dim; jmm++) {
        a = a - dot<scalar>(DphiiTinv.column(jmm),v)*aDphilmujT(n,i,jmm,jm,mujnj,L,Gt,rhoj0,r2,sw2);

    }

    if (i == n && j == jm) {
        a = a - dtDphi(n,DphiiTinv.T()*v,Gt,rhoj0,L,r2,sw2);
    }

    return a;
}

