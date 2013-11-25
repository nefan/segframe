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

% sample image and derivatives
function [sI,sD1,sD2] = linSampleI(I,D1I,D2I,ps,order)
    dim = size(ps,1);
    L = size(ps,2);
    
    sI = interpn(I,ps(1,:),ps(2,:),'linear',0);

    sD1 = zeros(L,dim);
    for i=1:dim
        sD1(:,i) = interpn(D1I(:,:,i),ps(1,:),ps(2,:),'linear',0);
    end

    sD2 = [];    
    if order >= 1
        sD2 = zeros(L,dim,dim);
        for i=1:dim
            for j=1:dim
                sD2(:,i,j) = interpn(D2I(:,:,i,j),ps(1,:),ps(2,:),'linear',0);
            end
        end        
    end
end