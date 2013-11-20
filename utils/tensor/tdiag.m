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

function T = tdiag(T1,d1,d2)
%
% take diagonal of tensor T1 over dimensions d1 and d2, d1<d2 with result in d1
%

size1 = size(T1);

assert(d1<d2);
assert(size1(d1)==size1(d2));

s1 = [1:(d1-1) (d1+1):(d2-1) (d2+1):ndims(T1)];

T1s = tshift(T1,[s1 d1 d2]);

T1s = reshape(T1s,[],size1(d1),size1(d1));

T = zeros(size(T1s,1),size1(d1));
for i=1:size1(d1)
    T(:,i) = T1s(:,i,i);
end

T = reshape(T,[size1(s1) size1(d1)]);
T = tshift(T,[1:(d1-1) ndims(T) d1:(ndims(T)-1)]);
