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

function T = tprod(T1,T2)
%
% tensor product between tensors T1 and T2
%

T1T = reshape(T1.T,[],1);
T2T = reshape(T2.T,1,[]);

TT = T1T*T2T; % do product

T = tensor(TT,[T1.dims T2.dims]);

if isfield(T1,'indices') && isfield(T2,'indices')
    T.indices = [T1.indices T2.indices];
end