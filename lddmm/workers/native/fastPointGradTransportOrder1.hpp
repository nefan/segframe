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

#ifndef FASTPOINTGRADTRANSPORTORDER1_H
#define FASTPOINTGRADTRANSPORTORDER1_H

#include "common.h"
#include "sla.h"
#include "gaussianKernels.hpp"
#include "order1.hpp"

const int dim = 3;
const int phioffset = 0;
const int Dphioffset = phioffset+dim;
const int muoffset = Dphioffset+dim*dim;
const int mujoffset = muoffset+dim;
int CSP = dim*(1+dim+1);
int gradCSP = 2*dim+2*dim*dim;

#define INDPHI(particle) gradCSP*particle+phioffset
#define INDDPHI(particle) gradCSP*particle+Dphioffset
#define INDMU(particle) gradCSP*particle+muoffset
#define INDMUJ(particle) gradCSP*particle+mujoffset

#endif // FASTPOINTGRADTRANSPORTORDER1_H
