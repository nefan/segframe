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

function q = getq(Gt,t,lddmmoptions)
%
% extract point positionis and derivatives
%

order = lddmmoptions.order;
dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim; % computations performed in cdim
dimq = lddmmoptions.dimq;
L = lddmmoptions.L;
R = lddmmoptions.R;
assert(R == 1);

x = reshape(deval(Gt,t),[],L);

if cdim == dim
    q = x(1:dimq,:);
else
    assert(dim == 2 && cdim == 3); % shift from 3d to 2d
    q = zeros(dimq,L);
    q(1:dim,:) = x(1:dim,:);    
    if order >= 1
        q(dim+(1:dim^2),:) = reshape(T23dTo2d(x(cdim+(1:cdim^2),:),lddmmoptions),dim^2,L); % q1
    end
    if order >= 2
        q(dim+dim^2+(1:dim^3),:) = reshape(T33dTo2d(x(cdim+cdim^2+(1:cdim^3),:),lddmmoptions),dim^3,L); % q2
    end
end

end