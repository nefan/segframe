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

function pointPath = getPointPathOrder1(moving,lddmmoptions)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim;
L = lddmmoptions.L;
R = lddmmoptions.R;
cCSP = lddmmoptions.cCSP;
% scales = lddmmoptions.scales;
% scaleweight = lddmmoptions.scaleweight;
% epsilon = lddmmoptions.epsilon;

ks = dkernelsGaussian(cdim);

    function Gt = pointPathOrder1(x)
        x = reshape(x,(1+dim)*dim,L);
        rhoj0 = x(dim+(1:dim^2),:);
        if dim ~= cdim
            rhoj0 = reshape(rhoj2dTo3dOrder1(rhoj0,lddmmoptions),cdim^2,L);
        end                

        initial = zeros(cCSP,L);
        initial(1:dim,:) = moving; % position
        initial(cdim+(1:cdim^2),:) = repmat(reshape(eye(cdim),cdim^2,1),1,L); % Dphi
        initial(cdim+cdim^2+(1:dim),:) = x(1:R*dim,:); % mu/rho0
        initial(cdim+cdim^2+cdim+(1:cdim^2),:) = -rhoj0;
        
        options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        [G] = getEvolutionEqs(lddmmoptions,rhoj0);
        Gt = ode45(G,[0 1],reshape(initial,L*cCSP,1),options); % matrix version
        assert(Gt.x(end) == 1); % if not, integration failed

    end

pointPath = @pointPathOrder1;

end