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

function gradU = getImageUL2(IM,IF,moving,imageoptions,dim,L)

scale = imageoptions.scale;
range = imageoptions.range;
order = imageoptions.order;
h = imageoptions.h;
switch order
    case 0
        dimX1 = dim;
    case 1
        dimX1 = dim+dim^2;
    case 2
        dimX1 = 2*dim+dim^2;
    otherwise
        assert(false);
end
movingTransform = @(x) x;
fixedTransform = @(x) x;
if isfield(imageoptions,'movingTransform')
    movingTransform = imageoptions.movingTransform;
end
if isfield(imageoptions,'fixedTransform')
    fixedTransform = imageoptions.fixedTransform;
end

    function [IG,D1IG,D2IG] = smooth(I)
        % Gaussian filter
        assert(dim == 2); % only 2D filter implemented
        Gfilter = fspecial('gaussian',range,scale);

        IG = imfilter(I,Gfilter,'same'); % smooth
        D1IG = zeros([size(I) dim]); % differences
        D2IG = zeros([size(IF) dim dim]);
        for i=1:dim
            diffIGi = diff(IG,1,i);
            ii = zeros(1,dim); ii(i) = 1;
            D1IG(1+ii(1):end-ii(1),1+ii(2):end-ii(2),i) = 0.5*(diffIGi(1+ii(1):end,1+ii(2):end)+diffIGi(1:(end-ii(1)),1:(end-ii(2))));
            for j=1:dim
                ij = zeros(1,dim); ij(j) = 1;
                diffIGij = diff(D1IG(1+ii(1):end-ii(1),1+ii(2):end-ii(2),i),1,j);
                D2IG(1+ii(1)+ij(1):end-ii(1)-ij(1),1+ii(2)+ij(2):end-ii(2)-ij(2),i,j) = 0.5*(diffIGij(1+ij(1):end,1+ij(2):end)+diffIGij(1:(end-ij(1)),1:(end-ij(2))));    
            end
        end
    end

    function [sI,sD1,sD2] = sample(I,D1I,D2I,ps)
        sI = interpn(I,ps(1,:),ps(2,:));
        sD1 = zeros(L,dim);
        sD2 = zeros(L,dim,dim);
        for i=1:dim
            sD1(:,i) = interpn(D1I(:,:,i),ps(1,:),ps(2,:));
            for j=1:dim
                sD2(:,i,j) = interpn(D2I(:,:,i,j),ps(1,:),ps(2,:));
            end
        end        
    end

% smooth and sample IF
[IFG,D1IFG,D2IFG] = smooth(IF);
[IMG,D1IMG,D2IMG] = smooth(IM);
[sIFG,sD1IFG,sD2IFG] = sample(IFG,D1IFG,D2IFG,moving);

    function [y,v] = lgradU(x)
        %
        % gradient, image matching
        %                
   
        x = vtot(x,[dimX1 L]);
        q0_a_i = x(1:cdim,:);
        switch order
            case 1
                q1_a__bi = reshape(x((cdim+1):(cdim+cdim^2),:),cdim,cdim,L);
            otherwise
                assert(false);
        end
        y = 0;
        v = []; % zeros(dimX1,L); gradient not implemented yet

        ps = x(1:dim,:);
        [sIMG,sD1IMG,sD2IMG] = sample(IMG,D1IMG,D2IMG,q0_a_i);
        
        % measure similarity
        y = h*(sIFG-sIMG).^2+h^3/12*(sD1IFG-tproddiag(sD1IMG,1,q1_a__bi,1)).^2;
        if order > 1
            assert(false); % add second order terms
        end
        
%         y = y/L;
%         v = reshape(v,dimX1*L,1)/L;
    end

gradU = @lgradU;

end