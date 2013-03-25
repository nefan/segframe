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

function [grid1 gridt] = shootgridDKernels(x,grid0,Gtt,lddmmoptions,varargin)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim;
L = lddmmoptions.L;
R = lddmmoptions.R;
CSP = lddmmoptions.CSP;
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;
Ngrid = numel(grid0{1});

x = reshape(x,(1+dim)*dim,L);
rhoj0 = reshape(x(dim+(1:dim^2),:),dim^2*L,1);
if dim ~= cdim
    rhoj0 = rhoj2dTo3dOrder1(rhoj0,lddmmoptions);
end  

backwards = false;
if size(varargin,2) > 0
    backwards = varargin{1};
end

function M = Ks(x,y,r,sw)
    M = 1/sw^2*exp(-sum((x-y).^2)/r^2);
end
function M = N1Ks(x,y,r,sw)
    M = -2/(sw^2*r^2)*exp(-sum((x-y).^2)/r^2)*(x-y);
end

function dgrid = Gc(t,gridt) % wrapper for C version
    t = intTime(t,backwards,lddmmoptions);
    
    Gt = deval(Gtt,t);
           
    dgrid = fastPointTransportOrder1(t,gridt,Gt,rhoj0,L,R,cdim,scales.^2,scaleweight.^2);

    dgrid = intResult(dgrid,backwards,lddmmoptions);
    
    % debug
    global testC
    if testC
        dgrid2 = G(tt,gridt);
        assert(norm(dgrid-dgrid2) < 10e-12);
    end
end

function dgrid = G(t,grid) % slooow version
    t = intTime(t,backwards,lddmmoptions);    
    
    grid = reshape(grid,3,Ngrid);
    
    dgrid = zeros(size(grid));
    
    Gt = reshape(deval(Gtt,t),CSP,L);
    
    for n = 1:Ngrid % grid point        
        for i = 1:L % particle
            xi = Gt(1:dim,i);
            Dphii = reshape(Gt(dim+(1:dim^2),i),dim,dim);            
            
            for si = 1:R % scale
                mui = Gt(dim+dim^2+(1:dim),i);
                
                % position
                dgrid(1:dim,n) = dgrid(1:dim,n)...
                    +Ks(xi,grid(1:dim,n),scales(si),scaleweight(si))*mui;
                for j = 1:dim
                    muij = Dphii'\rhoj0(dim*(j-1)+(1:dim),i);
                    
                    dgrid(1:dim,n) = dgrid(1:dim,n)...
                        +N1Ks(xi,grid(1:dim,n),scales(si),scaleweight(si))'*Dphii(:,j)*muij;
                end
            end
        end
    end 
    
    dgrid = intResult(dgrid,backwards,lddmmoptions);
    
    dgrid = reshape(dgrid,3*Ngrid,1);    
end

siz = size(grid0{1});
g0 = reshape([reshape(grid0{1},1,Ngrid); reshape(grid0{2},1,Ngrid); reshape(grid0{3},1,Ngrid)],3*Ngrid,1);
gridt = ode45(@Gc,[0 1],g0);
% gridt = ode45(@G,[0 1],g0);
g1 = reshape(deval(gridt,1),3,Ngrid);
grid1 = cell(1);
grid1{1} = reshape(g1(1,:),siz);
grid1{2} = reshape(g1(2,:),siz);
grid1{3} = reshape(g1(3,:),siz);

end