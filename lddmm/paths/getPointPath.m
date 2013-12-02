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

function pointPath = getPointPath(lddmmoptions)

order = lddmmoptions.order;

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim; % computations performed in cdim
L = lddmmoptions.L;
R = lddmmoptions.R;
CSP = lddmmoptions.CSP;
cCSP = lddmmoptions.cCSP;

function Gt = lpointPath(x, varargin)   
    tend = 1;
    if size(varargin,2) > 0
        tend = varargin{1};
    end
    x = reshape(x,CSP,L);
    

    if dim == cdim
        cx = x;        
    else
        assert(dim == 2 && cdim == 3); % shift from 2d to 3d
        cx = zeros(cCSP,L);
        cx(1:dim,:) = x(1:dim,:); % q0
        
        switch order
            case 0
                cx(cdim+(1:cdim:cdim*R),:) = x(dim+(1:dim:dim*R),:);
                cx(cdim+(2:cdim:cdim*R),:) = x(dim+(2:dim:dim*R),:);
            case 1
                cx(cdim+(1:cdim^2),:) = reshape(T22dTo3d(x(dim+(1:dim^2),:),lddmmoptions),cdim^2,L); % q1
                cx(cdim+cdim^2,:) = 1; % add one in last dim
                cx(cdim+cdim^2+(1:dim),:) = x(dim+dim^2+(1:dim),:); % mu0
                cx(2*cdim+cdim^2+(1:cdim^2),:) = -reshape(T22dTo3d(x(2*dim+dim^2+(1:dim^2),:),lddmmoptions),cdim^2,L); % mu1 !!!! note the minus
            case 2
                cx(cdim+(1:cdim^2),:) = reshape(T22dTo3d(x(dim+(1:dim^2),:),lddmmoptions),cdim^2,L); % q1
                cx(cdim+cdim^2,:) = 1; % add one in last dim
                cx(cdim+cdim^2+(1:cdim^3),:) = reshape(T32dTo3d(x(dim+dim^2+(1:dim^3),:),lddmmoptions),cdim^3,L); % q2
                cx(cdim+cdim^2+cdim^3+(1:dim),:) = x(dim+dim^2+dim^3+(1:dim),:); % mu0
                cx(2*cdim+cdim^2+cdim^3+(1:cdim^2),:) = reshape(T22dTo3d(x(2*dim+dim^2+dim^3+(1:dim^2),:),lddmmoptions),cdim^2,L); % mu1                
                cx(2*cdim+2*cdim^2+cdim^3+(1:cdim^3),:) = reshape(T32dTo3d(x(2*dim+2*dim^2+dim^3+(1:dim^3),:),lddmmoptions),cdim^3,L); % mu2
        end
    end
    
    options = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [G] = getEvolutionEqs(lddmmoptions);
    Gt = ode45(G,[0 tend],reshape(cx,L*cCSP,1),options); % matrix version
    assert(Gt.x(end) == tend); % if not, integration failed    
end

pointPath = @lpointPath;

end