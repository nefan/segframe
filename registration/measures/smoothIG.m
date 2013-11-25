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

%  smooth with Gaussian filter and take derivatives
function [IG,D1IG,D2IG] = smoothIG(I,imageoptions)
    scale = imageoptions.scale;
    range = imageoptions.range;
    order = imageoptions.order;
    dim = imageoptions.dim;
    
    assert(ndims(I) == 2); % only 2D filter implemented
    Gfilter = fspecial('gaussian',range,scale);

    IG = imfilter(I,Gfilter,'same'); % smooth

    D1IG = zeros([size(I) dim]); % differences
    for i=1:dim
        diffIGi = diff(IG,1,i);
        ii = zeros(1,dim); ii(i) = 1;
        D1IG(1+ii(1):end-ii(1),1+ii(2):end-ii(2),i) = 0.5*(diffIGi(1+ii(1):end,1+ii(2):end)+diffIGi(1:(end-ii(1)),1:(end-ii(2))));
    end

    D2IG = [];
    if order >= 1
        D2IG = zeros([size(IG) dim dim]);
        for i=1:dim
            ii = zeros(1,dim); ii(i) = 1;
            for j=1:dim
                ij = zeros(1,dim); ij(j) = 1;
                diffIGij = diff(D1IG(1+ii(1):end-ii(1),1+ii(2):end-ii(2),i),1,j);
                D2IG(1+ii(1)+ij(1):end-ii(1)-ij(1),1+ii(2)+ij(2):end-ii(2)-ij(2),i,j) = 0.5*(diffIGij(1+ij(1):end,1+ij(2):end)+diffIGij(1:(end-ij(1)),1:(end-ij(2))));    
            end
        end
    end
end