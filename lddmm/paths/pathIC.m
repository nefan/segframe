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

function y = pathIC(x,moving,lddmmoptions)

order = lddmmoptions.order;

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim; % computations performed in cdim
L = lddmmoptions.L;
R = lddmmoptions.R;
assert(R == 1);

x = reshape(x,[],L);
switch order
    case 0
        % variables: mu0
        q0 = moving;
        y = [q0; x];
    case 1
        q0 = moving;
        q1 = repmat(reshape(eye(dim),dim^2,1),1,L);
        y = [q0; q1; x];
    case 2
        q0 = moving;
        q1 = repmat(reshape(eye(dim),dim^2,1),1,L);
        q2 = zeros(dim^3,L);
        y = [q0; q1; q2; x];                
end

end