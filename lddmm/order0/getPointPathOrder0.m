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

function pointPath = getPointPathOrder0(moving,lddmmoptions)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim; % computations performed in cdim
L = lddmmoptions.L;
R = lddmmoptions.R;
cCSP = lddmmoptions.cCSP;

    function rhot = pointPathOrder0(x)
        if dim == cdim
            rho0 = [moving; reshape(x,dim*R,L)];
        else
            assert(dim == 2 && cdim == 3); % shift from 2d to 3d
            rho0 = zeros(cdim+R*cdim,L);
            rho0(1:dim,:) = moving;
            xx = reshape(x,dim*R,L);
            rho0(cdim+(1:cdim:cdim*R),:) = xx(1:dim:dim*R,:);
            rho0(cdim+(2:cdim:cdim*R),:) = xx(2:dim:dim*R,:);
        end
        
        options = odeset('RelTol',1e-6,'AbsTol',1e-6);   
        [G] = getEvolutionEqs(lddmmoptions);
        rhot = ode45(G,[0 1],reshape(rho0,L*cCSP,1),options);
        assert(rhot.x(end) == 1); % if not, integration failed     
    end

pointPath = @pointPathOrder0;

end
