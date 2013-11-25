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

% smooth and sample IF
[IFG,D1IFG,D2IFG] = smoothIG(IF,imageoptions);
[IMG,D1IMG,D2IMG] = smoothIG(IM,imageoptions);
[sIFG,sD1IFG,sD2IFG] = linSampleI(IFG,D1IFG,D2IFG,fixedTransform(moving),order);

    function [y,v] = lgradU(x)
        %
        % gradient, image matching
        %                
   
        x = vtot(x,[dimX1 L]);
        q0_a_i = x(1:dim,:);
        
        % sample
        [sIMG,sD1IMG,sD2IMG] = linSampleI(IMG,D1IMG,D2IMG,movingTransform(q0_a_i),order);

        if order == 1
            q1_a__bi = reshape(x((dim+1):(dim+dim^2),:),dim,dim,L);
        end
        if order == 2
            assert(false);
        end
                
        % measure similarity
        ys0 = (sIFG-sIMG)';
        y = h*sum(ttov(ys0).^2);
        if order == 1
            ys1 = sD1IFG-tcntr(tproddiag(sD1IMG,1,q1_a__bi,3),2,3);
            y = y + h^3/12*sum(ttov(ys1).^2);
        end
        if order == 2
            assert(false); % add second order terms
        end

        % gradient
        v = zeros(dimX1,L);
        v(1:dim,:) = tshift(-2*h*tproddiag(ys0,1,sD1IMG,1),[2 1]);
        
        if order == 1
            v(1:dim,:) = v(1:dim,:)+tshift(-2*h*tcntr(tproddiag(ys1,1,tcntr(tproddiag(sD2IMG,1,q1_a__bi,3),2,5),1),2,4),[2 1]);
            v(dim+(1:dim^2),:) = reshape(tshift(-2*h^3/12*tdiag(tproddiag(ys1,1,tdiag(tproddiag(sD1IMG,1,q1_a__bi,3),2,3),1),2,4),[3 2 1]),dim^2,L);
        end
        if order == 2
            assert(false); % add terms
        end        
               
        y = y/L;
        v = reshape(v,dimX1*L,1)/L;
    end

gradU = @lgradU;

end