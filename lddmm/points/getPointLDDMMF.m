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

function gradF = getPointLDDMMF(gradU,moving,methods,lddmmoptions)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim;
L = lddmmoptions.L;
R = lddmmoptions.R;
CSP = lddmmoptions.CSP;
dimq = lddmmoptions.dimq;
order = lddmmoptions.order;
energyweight = lddmmoptions.energyweight;

pointPath = methods.pointPath;
gradTransport = methods.gradTransport;
pathEnergy = methods.pathEnergy;

    function [y v] = lgradF(x)
        %
        % gradient by pull back
        %       
        
        % forward shot and path energy
        Gt = pointPath(pathIC(x,moving,lddmmoptions)); % Gt is a solver structure
        shot = getq(Gt,1,lddmmoptions);

        % warning: order1 does not return gradient for x-parts
        % warning: order0 adds energy part during backwards gradient transport
        if order == 0        
            Epath = pathEnergy(x,Gt); 
        else
            [Epath vE] = pathEnergy(x,Gt); 
        end
                
        % gradient at endpoint
        [U,v1] = gradU(reshape(shot,dimq*L,1));        
        y = energyweight*[Epath; U];    
        if nargout == 1 || numel(v1) == 0 % gradient not needed or not provided by gradU
            return
        end        
        v1 = reshape(v1,[],L);
        w1 = reshape([v1; zeros(lddmmoptions.CSP-size(v1,1),L)],CSP*L,1);
        
        % pull back
        % unfortunately, a bit unclean due to the incorporation of the
        % energy term in the order0 transport
        if order == 0
            vU = gradTransport(energyweight(2)*w1,x,Gt);
            v = vU;
        else
            assert(order <= 1);            
            vU = gradTransport(w1,x,Gt);
            v = energyweight(1)*vE+energyweight(2)*vU;
        end  
        
        if isfield(methods,'prior')
            [py,pv] = methods.prior(x);
            y = y + py;
            v = v + pv;
        end
        
        % no point movement at this point
        v = reshape(v,CSP,L);
        v = v(dimq+(1:CSP-dimq),:);
        v = reshape(v,(CSP-dimq)*L,1);
    end

gradF = @lgradF;

end