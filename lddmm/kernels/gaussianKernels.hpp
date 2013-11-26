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

#ifndef GAUSSIANKERNELS_H
#define GAUSSIANKERNELS_H

using namespace sla;

namespace gs {

inline scalar gamma(const Vector3<scalar> &xmy, const scalar r2, const scalar sw2) {
    return 1.0/sw2*exp(-dot<scalar>(xmy,xmy)/r2);
}
inline scalar d1gamma(const Vector3<scalar> &xmy, const scalar r2, const scalar sw2) {
    return -1.0/(sw2*r2)*exp(-dot<scalar>(xmy,xmy)/r2);
}
inline scalar d2gamma(const Vector3<scalar> &xmy, const scalar r2, const scalar sw2) {
    return 1.0/(sw2*r2*r2)*exp(-dot<scalar>(xmy,xmy)/r2);
}
inline scalar d3gamma(const Vector3<scalar> &xmy, const scalar r2, const scalar sw2) {
    return -1.0/(sw2*r2*r2*r2)*exp(-dot<scalar>(xmy,xmy)/r2);
}
inline scalar d4gamma(const Vector3<scalar> &xmy, const scalar r2, const scalar sw2) {
    const scalar r4 = r2*r2;
    return 1.0/(sw2*r4*r4)*exp(-dot<scalar>(xmy,xmy)/r2);
}

// Hermite function
inline scalar He(const int n, const scalar x) {
    scalar v;

    switch(n) {
        case 0:
            v = 1;
            break;
        case 1:
            v = x;
            break;
        case 2:
            v = x*x-1;
            break;
        case 3:
            v = pow(x,3)-3*x;
            break;
        case 4:
            v = pow(x,4)-6*x*x+3;
            break;
        case 5:
            v = pow(x,5)-10*pow(x,3)+15*x;
            break;
        default:
            assert(false);
    }

    return v;
}

// Gaussian derivative
inline scalar DaKs(Vector3<int> &da, const Vector3<scalar> &xmy, const scalar r, const scalar ks) {
    const Vector3<scalar> z = sqrt(2)/r*xmy;

    return pow(-sqrt(2)/r,sum(da))*He(da[0],z[0])*He(da[1],z[1])*He(da[2],z[2])*ks;
}

};

#endif // GAUSSIANKERNELS_H
