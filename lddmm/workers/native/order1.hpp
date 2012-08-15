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

#ifndef ORDER1_H
#define ORDER1_H

#include "sla.h"
#include "gaussianKernels.hpp"

using namespace sla;
using namespace gs;

namespace kd {

inline scalar Ks(const Vector3<scalar> &xmy, const scalar r2, const scalar sw2) {
    return gamma(xmy,r2,sw2);
}

inline Vector3<scalar> N1Ks(const Vector3<scalar> &xmy, const scalar r2, const scalar sw2) {
    return 2.0*d1gamma(xmy,r2,sw2)*xmy;
}

inline Vector3<scalar> N2Ks(const Vector3<scalar> &xmy, const scalar r2, const scalar sw2) {
    return -2.0*d1gamma(xmy,r2,sw2)*xmy;
}

inline Vector3<scalar> D1N2Ks(const Vector3<scalar> &xmy, const Vector3<scalar> &v, const scalar r2, const scalar sw2) {
    return -2.0*d1gamma(xmy,r2,sw2)*v
        -4.0*d2gamma(xmy,r2,sw2)*dot<scalar>(xmy,v)*xmy;
}
inline Vector3<scalar> D1N2KsT(const Vector3<scalar> &xmy, const Vector3<scalar> &v, const scalar r2, const scalar sw2) {
    return D1N2Ks(xmy,v,r2,sw2);
}

inline Vector3<scalar> D2N2Ks(const Vector3<scalar> &xmy, const Vector3<scalar> &v, const scalar r2, const scalar sw2) {
    return 2.0*d1gamma(xmy,r2,sw2)*v
        +4.0*d2gamma(xmy,r2,sw2)*dot<scalar>(xmy,v)*xmy;
}
inline Vector3<scalar> D2N2KsT(const Vector3<scalar> &xmy, const Vector3<scalar> &v, const scalar r2, const scalar sw2) {
    return D2N2Ks(xmy,v,r2,sw2);
}

inline Vector3<scalar> D2D1N2Ks(const Vector3<scalar> &xmy, const Vector3<scalar> &a, const Vector3<scalar> &v, const scalar r2, const scalar sw2) {
    return 4.0*d2gamma(xmy,r2,sw2)*dot<scalar>(a,v)*xmy
        +4.0*d2gamma(xmy,r2,sw2)*dot<scalar>(xmy,v)*a
        +4.0*d2gamma(xmy,r2,sw2)*dot<scalar>(xmy,a)*v
        +8.0*d3gamma(xmy,r2,sw2)*dot<scalar>(xmy,a)*dot<scalar>(xmy,v)*xmy;
}
inline Vector3<scalar> D2D1N2KsT(const Vector3<scalar> &xmy, const Vector3<scalar> &a, const Vector3<scalar> &v, const scalar r2, const scalar sw2) {
    return D2D1N2Ks(xmy,a,v,r2,sw2);
}

};

#endif // ORDER1_H

