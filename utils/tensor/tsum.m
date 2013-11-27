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

function T = tsum(T1,T2)
%
% elementwise sum of tensors
% indices are removed if they do not correspond
%

assert(isequal(T1.dims,T2.dims));

T = tensor(T1.T+T2.T,T1.dims);

if isfield(T1,'indices') && isfield(T2,'indices') ...
        && isequal(T1.indices,T2.indices)
    T.indices = T1.indices;
end