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

function gradU = getPointGradU(fixed,dim,L,order,varargin)

weights = [];
if size(varargin,2) > 0
    weights = varargin{1};
end

assert(size(fixed,2) == L);
if order == 0
    assert(size(fixed,1) == dim);
else
    assert(size(fixed,1) == dim+dim^2);
end

    function [y v] = lgradU(x)
        %
        % gradient, point match
        %        
        
        % initial value
        w = x-reshape(fixed,[],1);
        if ~isempty(weights)
            w = w.*repmat(weights,L,1);
        end
        y = sum(w.^2);
        v = reshape(2*w,[],1);

    end

gradU = @lgradU;

end        
