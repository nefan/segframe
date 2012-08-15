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

function M = dkernelDerivativeMatrix(rho0,rhoj0,lddmmoptions)

cdim = lddmmoptions.cdim;
L = lddmmoptions.L;
R = lddmmoptions.R;
CSP = lddmmoptions.CSP;
gradCSP = 2*cdim+2*cdim^2;
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;
epsilon = lddmmoptions.epsilon;

assert(R == 1)
s = 1; % only one scale

ks = dkernelsGaussian(cdim);

function xi = x(i,Gt)
    xi = Gt(1:cdim,i);
end
function mui = mu(i,Gt)
    mui = Gt(cdim+cdim^2+(1:cdim),i);
end
function muji = muj(i,j,Gt)
    muji = Dphi(i,Gt)'\rhoj0(cdim*(j-1)+(1:cdim),i);
end
function Dphi = Dphi(i,Gt)
    Dphi = reshape(Gt(cdim+(1:cdim*cdim),i),cdim,cdim);
end
function Dphil = Dphil(i,l,Gt)
    Dphil = Gt(cdim+cdim*(l-1)+(1:cdim),i);
end
function dtDphi = dtDphi(n,Gt)
    dtDphi = zeros(cdim,cdim);
    xn = x(n,Gt);

    for c = 1:cdim
        for i = 1:L
            xi = x(i,Gt);
            dtDphi(:,c) = dtDphi(:,c)...
                +ks.N2Ks(xi,xn,scales(s),scaleweight(s))'*Dphil(n,c,Gt)*mu(i,Gt);
            for j = 1:cdim
                dtDphi(:,c) = dtDphi(:,c)...
                    +(ks.D1N2Ks(xi,xn,scales(s),scaleweight(s))*Dphil(i,j,Gt))'*Dphil(n,c,Gt)*muj(i,j,Gt);
            end        
        end
    end
end

function a = aphiphi(n,i,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = 2*ks.d1gamma(xi,xn,scales(s),scaleweight(s))*mu(i,Gt)*(xi-xn)';
    
    for jm = 1:cdim
        a = a + muj(i,jm,Gt)*Dphil(i,jm,Gt)'...
            *(4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*(xi-xn)'...
             +2*ks.d1gamma(xi,xn,scales(s),scaleweight(s))*eye(cdim));
    end
    
    if i == n
        for im = 1:L
            xim = x(im,Gt);
            a = a ...
                -2*ks.d1gamma(xim,xn,scales(s),scaleweight(s))*mu(im,Gt)*(xim-xn)';
            
            for jm = 1:cdim
                a = a - muj(im,jm,Gt)*Dphil(im,jm,Gt)'...
                        *(4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)*(xim-xn)'...
                         +2*ks.d1gamma(xim,xn,scales(s),scaleweight(s))*eye(cdim));
            end
        end
    end
end
function a = aphiDphi(n,i,l,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = muj(i,l,Gt)*ks.N1Ks(xi,xn,scales(s),scaleweight(s))';
end
function a = aphimu(n,i,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = ks.gamma(xi,xn,scales(s),scaleweight(s));
end
function a = aphimuj(n,i,j,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = ks.N1Ks(xi,xn,scales(s),scaleweight(s))'*Dphil(i,j,Gt);
end
function a = aDphiphi(n,i,l,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = -mu(i,Gt)*Dphil(n,l,Gt)'...
         *(4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*(xi-xn)'...
          +2*ks.d1gamma(xi,xn,scales(s),scaleweight(s))*eye(cdim));
    
    for jm = 1:cdim
        a = a...
            -muj(i,jm,Gt)*Dphil(n,l,Gt)'...
             *(4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*Dphil(i,jm,Gt)*(xi-xn)'...
              +8*ks.d3gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*(xi-xn)'*Dphil(i,jm,Gt)*(xi-xn)'...
              +4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)'*Dphil(i,jm,Gt)*eye(cdim)...
              +4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*Dphil(i,jm,Gt)');
    end
    
    if i == n
        for im = 1:L
            xim = x(im,Gt);
            a = a...
                +mu(im,Gt)*Dphil(n,l,Gt)'...
                 *(4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)*(xim-xn)'...
                  +2*ks.d1gamma(xim,xn,scales(s),scaleweight(s))*eye(cdim));
            
            for jm = 1:cdim
                a = a...
                    +muj(im,jm,Gt)*Dphil(n,l,Gt)'...
                     *(4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*Dphil(im,jm,Gt)*(xim-xn)'...
                      +8*ks.d3gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)*(xim-xn)'*Dphil(im,jm,Gt)*(xim-xn)'...
                      +4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)'*Dphil(im,jm,Gt)*eye(cdim)...
                      +4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)*Dphil(im,jm,Gt)');
            end
        end
    end
end
function a = aDphiDphi(n,i,l,k,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = muj(i,k,Gt)*Dphil(n,l,Gt)'*ks.D1N2Ks(xi,xn,scales(s),scaleweight(s));
    
    if i == n && l == k
        for im = 1:L
            xim = x(im,Gt);
            a = a + mu(im,Gt)*ks.N2Ks(xim,xn,scales(s),scaleweight(s))';
            
            for jm = 1:cdim
                a = a + muj(im,jm,Gt)*(ks.D1N2Ks(xim,xn,scales(s),scaleweight(s))*Dphil(im,jm,Gt))';
            end
        end
    end
end
function a = aDphimu(n,i,l,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = ks.N2Ks(xi,xn,scales(s),scaleweight(s))'*Dphil(n,l,Gt);
end
function a = aDphimuj(n,i,l,j,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = (ks.D1N2Ks(xi,xn,scales(s),scaleweight(s))*Dphil(i,j,Gt))'*Dphil(n,l,Gt);
end
function a = amuphi(n,i,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = mu(n,Gt)'*mu(i,Gt)*(4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*(xi-xn)' ...
        +2*ks.d1gamma(xi,xn,scales(s),scaleweight(s))*eye(cdim));
    
    for jm = 1:cdim
        a = a ...
            -(muj(n,jm,Gt)'*mu(i,Gt)-mu(n,Gt)'*muj(i,jm,Gt))...
             *(4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*Dphil(n,jm,Gt)*(xi-xn)'...
              +8*ks.d3gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*(xi-xn)'*Dphil(n,jm,Gt)*(xi-xn)'...
              +4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)'*Dphil(n,jm,Gt)*eye(cdim)...
              +4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*Dphil(n,jm,Gt)');
          
        for j = 1:cdim
            a = a ...
                -(muj(n,jm,Gt)'*muj(i,j,Gt))...
                 *(8*ks.d3gamma(xi,xn,scales(s),scaleweight(s))*Dphil(i,j,Gt)*(xi-xn)'*Dphil(n,jm,Gt)*(xi-xn)'...
                  +4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*Dphil(i,j,Gt)*Dphil(n,jm,Gt)'...
                  +8*ks.d3gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*Dphil(i,j,Gt)'*Dphil(n,jm,Gt)*(xi-xn)'...
                  +4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*Dphil(i,j,Gt)'*Dphil(n,jm,Gt)*eye(cdim)...
                  +8*ks.d3gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)'*Dphil(i,j,Gt)*Dphil(n,jm,Gt)*(xi-xn)'...
                  +4*ks.d2gamma(xi,xn,scales(s),scaleweight(s))*Dphil(n,jm,Gt)*Dphil(i,j,Gt)'...
                  +16*ks.d4gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)'*Dphil(i,j,Gt)*(xi-xn)*(xi-xn)'*Dphil(n,jm,Gt)*(xi-xn)'...
                  +8*ks.d3gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*(xi-xn)'*Dphil(n,jm,Gt)*Dphil(i,j,Gt)'...
                  +8*ks.d3gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)'*Dphil(i,j,Gt)*(xi-xn)'*Dphil(n,jm,Gt)*eye(cdim)...
                  +8*ks.d3gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)'*Dphil(i,j,Gt)*(xi-xn)*Dphil(n,jm,Gt)');
        end
    end
    
    if i == n
        for im = 1:L
            xim = x(im,Gt);
            a = a ...
                -mu(n,Gt)'*mu(im,Gt)*(4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)*(xim-xn)' ...
                +2*ks.d1gamma(xim,xn,scales(s),scaleweight(s))*eye(cdim));
            
                for jm = 1:cdim
                    a = a...
                        +(muj(n,jm,Gt)'*mu(im,Gt)-mu(n,Gt)'*muj(im,jm,Gt))...
                         *(4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*Dphil(n,jm,Gt)*(xim-xn)'...
                          +8*ks.d3gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)*(xim-xn)'*Dphil(n,jm,Gt)*(xim-xn)'...
                          +4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)'*Dphil(n,jm,Gt)*eye(cdim)...
                          +4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)*Dphil(n,jm,Gt)');

                    for j = 1:cdim
                        a = a...
                            +(muj(n,jm,Gt)'*muj(im,j,Gt))...
                             *(8*ks.d3gamma(xim,xn,scales(s),scaleweight(s))*Dphil(im,j,Gt)*(xim-xn)'*Dphil(n,jm,Gt)*(xim-xn)'...
                              +4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*Dphil(im,j,Gt)*Dphil(n,jm,Gt)'...
                              +8*ks.d3gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)*Dphil(im,j,Gt)'*Dphil(n,jm,Gt)*(xim-xn)'...
                              +4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*Dphil(im,j,Gt)'*Dphil(n,jm,Gt)*eye(cdim)...
                              +8*ks.d3gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)'*Dphil(im,j,Gt)*Dphil(n,jm,Gt)*(xim-xn)'...
                              +4*ks.d2gamma(xim,xn,scales(s),scaleweight(s))*Dphil(n,jm,Gt)*Dphil(im,j,Gt)'...
                              +16*ks.d4gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)'*Dphil(im,j,Gt)*(xim-xn)*(xim-xn)'*Dphil(n,jm,Gt)*(xim-xn)'...
                              +8*ks.d3gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)*(xim-xn)'*Dphil(n,jm,Gt)*Dphil(im,j,Gt)'...
                              +8*ks.d3gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)'*Dphil(im,j,Gt)*(xim-xn)'*Dphil(n,jm,Gt)*eye(cdim)...
                              +8*ks.d3gamma(xim,xn,scales(s),scaleweight(s))*(xim-xn)'*Dphil(im,j,Gt)*(xim-xn)*Dphil(n,jm,Gt)');
                    end
                end
        end
    end    
end
function a = amuDphi(n,i,l,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = zeros(cdim);
    
    for jm = 1:cdim
        a = a...
            -4*(muj(n,jm,Gt)'*muj(i,l,Gt))...
             *(ks.d2gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)'*Dphil(n,jm,Gt)*eye(cdim)...
              +ks.d2gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*Dphil(n,jm,Gt)'...
              +ks.d2gamma(xi,xn,scales(s),scaleweight(s))*Dphil(n,jm,Gt)*(xi-xn)'...
              +2*ks.d3gamma(xi,xn,scales(s),scaleweight(s))*(xi-xn)*(xi-xn)'*Dphil(n,jm,Gt)*(xi-xn)');        
    end
    
    if i == n
        for im = 1:L
            xim = x(im,Gt);
            a = a - (muj(n,l,Gt)'*mu(im,Gt)-mu(n,Gt)'*muj(im,l,Gt))*ks.D2N2Ks(xim,xn,scales(s),scaleweight(s));
            
            for jm = 1:cdim
                a = a - (muj(n,l,Gt)'*muj(im,jm,Gt))*ks.D2D1N2Ksa(Dphil(im,jm,Gt),xim,xn,scales(s),scaleweight(s));
            end
        end
    end
end
function a = amumu(n,i,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = -ks.N2Ks(xi,xn,scales(s),scaleweight(s))*mu(n,Gt)';
    
    for jm = 1:cdim
        a = a - ks.D2N2Ks(xi,xn,scales(s),scaleweight(s))*Dphil(n,jm,Gt)*muj(n,jm,Gt)';
    end    
    
    if i == n
        for im = 1:L
            xim = x(im,Gt);
            a = a ...
                -ks.N2Ks(xim,xn,scales(s),scaleweight(s))*mu(im,Gt)';
            
            for jm = 1:cdim
                a = a + ks.D2N2Ks(xim,xn,scales(s),scaleweight(s))*Dphil(n,jm,Gt)*muj(im,jm,Gt)';
            end
        end
    end
end
function a = amumuj(n,i,j,Gt)
    xi = x(i,Gt);
    xn = x(n,Gt);
    a = + ks.D2N2Ks(xi,xn,scales(s),scaleweight(s))*Dphil(n,j,Gt)*mu(n,Gt)';
    
    for jm = 1:cdim
        a = a - ks.D2D1N2Ksa(Dphil(i,j,Gt),xi,xn,scales(s),scaleweight(s))*Dphil(n,jm,Gt)*muj(n,jm,Gt)';
    end
    
    if i == n
        for im = 1:L
            xim = x(im,Gt);
            a = a - ks.D2N2Ks(xim,xn,scales(s),scaleweight(s))*Dphil(n,j,Gt)*mu(im,Gt)';
            
            for jm = 1:cdim
                a = a - ks.D2D1N2Ksa(Dphil(im,jm,Gt),xim,xn,scales(s),scaleweight(s))*Dphil(n,j,Gt)*muj(im,jm,Gt)';
            end
        end
    end
end
function a = amujphi(n,i,j,Gt)
    a = zeros(cdim);
    
    for jm = 1:cdim
        ejm = zeros(cdim,1); ejm(jm) = 1;
        a = a - (Dphi(i,Gt)'\ejm)*muj(n,j,Gt)'*aDphiphi(n,i,jm,Gt);
    end
end
function a = amujDphi(n,i,j,l,Gt)
    el = zeros(cdim,1); el(l) = 1;
    a = zeros(cdim);
    
    for jm = 1:cdim
        ejm = zeros(cdim,1); ejm(jm) = 1;
        a = a - (Dphi(i,Gt)'\ejm)*muj(n,j,Gt)'*aDphiDphi(n,i,jm,l,Gt);        
    end
    
    if i == n                
        a = a + (Dphi(n,Gt)'\el)*(Dphi(n,Gt)'\(dtDphi(n,Gt)'*muj(n,j,Gt)))';            
    end
end
function a = amujmu(n,i,j,Gt)
    a = zeros(cdim);
    
    for jm = 1:cdim
        ejm = zeros(cdim,1); ejm(jm) = 1;
        a = a - (Dphi(i,Gt)'\ejm)*muj(n,j,Gt)'*aDphimu(n,i,jm,Gt);
    end   
end
function a = amujmuj(n,i,j,jm,Gt)
    a = zeros(cdim);
    
    for jmm = 1:cdim
        ejmm = zeros(cdim,1); ejmm(jmm) = 1;
        a = a - (Dphi(i,Gt)'\ejmm)*muj(n,j,Gt)'*aDphimuj(n,i,jmm,jm,Gt);
    end   
    
    if i == n && j == jm
        a = a - (Dphi(n,Gt)')\(dtDphi(n,Gt)');
    end
end

M.x = @x;
M.mu = @mu;
M.muj = @muj;
M.Dphi = @Dphi;
M.Dphil = @Dphil;
M.dtDphi = @dtDphi;
M.aphiphi = @aphiphi;
M.aphiDphi = @aphiDphi;
M.aphimu = @aphimu;
M.aphimuj = @aphimuj;
M.aDphiphi = @aDphiphi;
M.aDphiDphi = @aDphiDphi;
M.aDphimu= @aDphimu;
M.aDphimuj = @aDphimuj;
M.amuphi = @amuphi;
M.amuDphi = @amuDphi;
M.amumu = @amumu;
M.amumuj = @amumuj;
M.amujphi = @amujphi;
M.amujDphi = @amujDphi;
M.amujmu = @amujmu;
M.amujmuj = @amujmuj;

end