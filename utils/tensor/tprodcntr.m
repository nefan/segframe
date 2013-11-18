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

function T = tcntr(T1,d1,T2,d2)
%
% contract tensor product of T1 and T2 over dimensions d1 and d2
%

size1 = size(T1);
size2 = size(T2);

s1 = [1:(d1-1) (d1+1):ndims(T1)];
s2 = [1:(d2-1) (d2+1):ndims(T2)];
T1s = tshift(T1,[s1 d1]);
T2s = tshift(T2,[d2 s2]);

T1s = reshape(T1s,[],size1(d1));
T2s = reshape(T2s,size2(d2),[]);

T = T1s*T2s; % do product and contraction

T = reshape(T,[size1(s1) size2(s2)]);