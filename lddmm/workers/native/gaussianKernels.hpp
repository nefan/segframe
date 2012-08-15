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

#include "fastPointGradTransportOrder1.hpp"

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

};

#endif // GAUSSIANKERNELS_H
