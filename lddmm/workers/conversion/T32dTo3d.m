%
%  segframe, Copyright (C) 2009-2012, Stefan Sommer (sommer@diku.dk)
%  https://github.com/nefan/segframe.git
% 
%  This file is part of segframe.
% 
%  segframe is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  segframe is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with segframe.  If not, see <http://www.gnu.org/licenses/>.
%  

function M3d = T22dTo3d(M2d,lddmmoptions)

%
% shift from 2d to 3d representation, M^3xL
%
dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim; % computations performed in cdim
L = lddmmoptions.L;
R = lddmmoptions.R;
assert(R == 1);

assert(dim == 2 && cdim == 3); % shift from 2d to 3d
M2d = reshape(M2d,dim,dim,dim,L);
M3d = zeros(cdim,cdim,cdim,L);
M3d(1:dim,1:dim,1:dim,:) = M2d;
M3d = reshape(M3d,cdim^3*L,1);