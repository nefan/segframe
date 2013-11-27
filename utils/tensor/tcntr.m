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

function T = tcntr(T1,d1,d2)
%
% contract tensor T1 over dimensions d1 and d2
%

assert(d1<d2);
assert(T1.dims(d1) == T1.dims(d2));

s1 = [1:(d1-1) (d1+1):(d2-1) (d2+1):tndims(T1)];

T1s = tshift(T1,[s1 d1 d2]);

T1sT = reshape(T1s.T,[],T1.dims(d1),T1.dims(d1));

TsT = zeros(size(T1sT,1),1);
for i=1:T1.dims(d1)
    TsT(:) = TsT(:) + T1sT(:,i,i);
end

T = tensor(TsT,T1.dims(s1));

if isfield(T1,'indices')
    T.indices = T1.indices(s1);
end