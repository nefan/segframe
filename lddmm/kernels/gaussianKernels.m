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

function [Ks D1Ks D2Ks] = gaussianKernels();
% scalar Gaussian kernel and derivative

% function M = lK(x,y,r,sw)
%     M = sw^(-2)*exp(-sum((x-y).^2)/r^2)*eye(d);
% end

function M = lKs(x,y,r,sw)
    M = sw^(-2)*exp(-sum((x-y).^2)/r^2);
end

% function M = lD1K(x,y,r,sw)
%     M = -2/(sw^2*r^2)*exp(-sum((x-y).^2)/r^2)*(x-y);
% end

function M = lD1Ks(x,y,r,sw)
    M = -1/(sw^2*r^2)*exp(-sum((x-y).^2)/r^2);
end

function M = lD2Ks(x,y,r,sw)
    M = 1/(sw^2*r^4)*exp(-sum((x-y).^2)/r^2);
end

% K = @lK;
Ks = @lKs;
% D1K = @lD1K;
D1Ks = @lD1Ks;
D2Ks = @lD2Ks;

end