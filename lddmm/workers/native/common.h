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

#include "mex.h"
#include "math.h"

typedef double scalar;

#define DEBUG

#define INDYX(particle,coord_particle,direction,coord_direction,scale_direction) CSD*direction+dim*scale_direction+coord_direction+CSDL*(CSP*particle+coord_particle)
#define INDYP(particle,coord_particle,scale_particle,direction,coord_direction,scale_direction) CSD*direction+dim*scale_direction+coord_direction+CSDL*(CSP*particle+dim+dim*scale_particle+coord_particle)

#define INDRHOX(particle,coord_particle) CSP*particle+coord_particle
#define INDRHOP(particle,coord_particle,scale_particle) CSP*particle+dim+dim*scale_particle+coord_particle

#define INDGX(particle,coord_particle) CSP*particle+coord_particle
#define INDGDPHI(particle,i,j) CSP*particle+dim+dim*j+i
#define INDGMU(particle,coord_particle) CSP*particle+dim*(1+dim)+coord_particle
//#define INDGMUJ(particle,i,j) CSP*particle+dim*(1+dim+1)+dim*j+i
#define INDRHOJ(particle,i,j) dim*dim*particle+dim*j+i

// old stuff, should be removed
#define KKs(kbase,sw2) kbase/sw2
#define D1Ks(kbase,r2,sw2) -kbase/(sw2*r2)
#define D2Ks(kbase,r2,sw2) kbase/(sw2*r2*r2)

#ifdef DEBUG
#define mexAssert(t) if (!(t)) mexErrMsgIdAndTxt("assert:assert",__STRING(t));
#else
#define mexAssert(t)
#endif
