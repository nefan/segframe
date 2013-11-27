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
sIFG = tensor(sIFG,[L],'i');
sD1IFG = tensor(sD1IFG,[L dim],'ia');
sD2IFG = tensor(sD2IFG,[L dim dim],'iab');

    function [y,v] = lgradU(x)
        %
        % gradient, image matching
        %                
   
        x = reshape(x,dimX1,L);
        q0 = tensor(x(1:dim,:),[dim L],'ai');
        
        % sample
        [sIMG,sD1IMG,sD2IMG] = linSampleI(IMG,D1IMG,D2IMG,movingTransform(q0.T),order);
        sIMG = tensor(sIMG,[L],'i');
        sD1IMG = tensor(sD1IMG,[L dim],'ia');
        sD2IMG = tensor(sD2IMG,[L dim dim],'iab');

        if order == 1
                q1 = tensor(x((dim+1):(dim+dim^2),:),[dim dim L],'abi');
        end
        if order == 2
            assert(false);
        end
                
        % measure similarity
        ys0 = tsub(sIFG,sIMG);
        y = h/L*sum(ttov(ys0).^2);
        if order == 1
            ys1 = tsub(sD1IFG,tcntr(tproddiag(sD1IMG,1,q1,3),2,3));
            ys1.indices = 'ia';
            y = y + h^3/(12*L)*sum(ttov(ys1).^2);
        end
        if order == 2
            assert(false); % add second order terms
        end

        % gradient
        v0 = tshift(tscalar(-2*h/L,tproddiag(ys0,1,sD1IMG,1)),[2 1]);
        
        if order == 1
            v0 = tsum(v0,tshift(tscalar(-2*h/L,tcntr(tproddiag(ys1,1,tcntr(tproddiag(sD2IMG,1,q1,3),2,5),1),2,4)),[2 1]));
            v1 = tshift(tscalar(-2*h^3/(12*L),tdiag(tproddiag(ys1,1,tdiag(tproddiag(sD1IMG,1,q1,3),2,3),1),2,4)),[3 2 1]);
        end
        if order == 2
            assert(false); % add terms
        end        
               
        v = reshape([v0.T; reshape(v1.T,dim^2,L)],dimX1*L,1);
    end

gradU = @lgradU;

end