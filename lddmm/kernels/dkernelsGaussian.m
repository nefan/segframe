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

function fs = dkernelsGaussian(dim)

function M = Ks(x,y,r,sw)
    M = 1/sw^2*exp(-sum((x-y).^2)/r^2);
end
function M = N1Ks(x,y,r,sw)
    M = -2/(sw^2*r^2)*exp(-sum((x-y).^2)/r^2)*(x-y);
end
function M = N2Ks(x,y,r,sw)
    M = 2/(sw^2*r^2)*exp(-sum((x-y).^2)/r^2)*(x-y);
end
function M = D1N1Ks(x,y,r,sw)
    M = -2/(sw^2*r^2)*exp(-sum((x-y).^2)/r^2)*eye(dim)...
        +4/(sw^2*r^4)*exp(-sum((x-y).^2)/r^2)*(x-y)*(x-y)';
end
function M = D2N1Ks(x,y,r,sw)
    M = 2/(sw^2*r^2)*exp(-sum((x-y).^2)/r^2)*eye(dim)...
        -4/(sw^2*r^4)*exp(-sum((x-y).^2)/r^2)*(x-y)*(x-y)';
end
function M = D1N2Ks(x,y,r,sw)
    M = 2/(sw^2*r^2)*exp(-sum((x-y).^2)/r^2)*eye(dim)...
        -4/(sw^2*r^4)*exp(-sum((x-y).^2)/r^2)*(x-y)*(x-y)';
end
function M = D2N2Ks(x,y,r,sw)
    M = -2/(sw^2*r^2)*exp(-sum((x-y).^2)/r^2)*eye(dim)...
        +4/(sw^2*r^4)*exp(-sum((x-y).^2)/r^2)*(x-y)*(x-y)';
end
function M = D2D1N2Ksa(a,x,y,r,sw)
    M = 4/(sw^2*r^4)*exp(-sum((x-y).^2)/r^2)*(a*(x-y)'+(x-y)*a'+(x-y)'*a*eye(dim))...
        -8/(sw^2*r^6)*exp(-sum((x-y).^2)/r^2)*(x-y)'*a*(x-y)*(x-y)';
end

function v = He(n,x)
    switch n
        case 0
            v = 1;
        case 1
            v = x;
        case 2
            v = x^2-1;
        case 3
            v = x^3-3*x;
        case 4
            v = x^4-6*x^2+3;
        case 5
            v = x^5-10*x^3+15*x;
        otherwise
            assert(false);
    end
end

function v = DaKs(da,x,y,r,sw)
    assert(dim == 3);
    z = sqrt(2)*(x-y)/r;
    v = 1/sw^2*(-sqrt(2)/r)^sum(da)*He(da(1),z(1))*He(da(2),z(2))*He(da(3),z(3))*exp(-sum(z.^2)/2);
end

function [D0Ks,D1Ks,D2Ks,D3Ks] = TKs(q0,scales,scaleweight)
    R = length(scales);
    assert(R == 1);
    sl = 1;
    L = q0.dims(2);
    cdim = q0.dims(1);

    function v = Kf(i,j)
        v = Ks(q0.T(:,i),q0.T(:,j),scales(sl),scaleweight(sl));
    end
    function v = D1Kf(i,j,b)
        da = zeros(cdim,1); da(b) = 1;            
        v = DaKs(da,q0.T(:,i),q0.T(:,j),scales(sl),scaleweight(sl));
    end
    function v = D2Kf(i,j,b,g)
        da = zeros(cdim,1); da(b) = da(b)+1; da(g) = da(g)+1;
        v = DaKs(da,q0.T(:,i),q0.T(:,j),scales(sl),scaleweight(sl));
    end        
    function v = D3Kf(i,j,b,g,d)
        da = zeros(cdim,1); da(b) = da(b)+1; da(g) = da(g)+1;  da(d) = da(d)+1;
        v = DaKs(da,q0.T(:,i),q0.T(:,j),scales(sl),scaleweight(sl));
    end

    % compute kernel and derivatives        

    D0Ks = tensor(mmap(@Kf,L,L),[L L],'ij');
    if nargout >= 2
        D1Ks = tensor(mmap(@D1Kf,L,L,cdim),[L L cdim],'ijb');
    end
    if nargout >=3
        D2Ks = tensor(mmap(@D2Kf,L,L,cdim,cdim),[L L cdim cdim],'ijbg');
    end
    if nargout >= 4
        D3Ks = tensor(mmap(@D3Kf,L,L,cdim,cdim,cdim),[L L cdim cdim cdim],'ijbgd');
    end        
end

function [D0Ks,D1Ks,D2Ks,D3Ks] = TKsC(q0_a_i,scales,scaleweight)
    R = length(scales);
    assert(R == 1);
    sl = 1;
    L = q0.dims(2);
    cdim = q0.dims(1);

    % compute kernel and derivatives
    switch nargout
        case 1
            [Ks__ij] = gaussianTKsC(q0.T,scales.^2,scaleweight.^2,int64(cdim),int64(L),int64(R));
            D0Ks = tensor(Ks__ij,[L L],'ij');
        case 2
            [Ks__ij,D1Ks__ijb] = gaussianTKsC(q0.T,scales.^2,scaleweight.^2,int64(cdim),int64(L),int64(R));
            D0Ks = tensor(Ks__ij,[L L],'ij');
            D1Ks = tensor(D1Ks__ijb,[L L cdim],'ijb');
        case 3
            [Ks__ij,D1Ks__ijb,D2Ks__ijbg] = gaussianTKsC(q0.T,scales.^2,scaleweight.^2,int64(cdim),int64(L),int64(R));
            D0Ks = tensor(Ks__ij,[L L],'ij');
            D1Ks = tensor(D1Ks__ijb,[L L cdim],'ijb');
            D2Ks = tensor(D2Ks__ijbg,[L L cdim cdim],'ijbg');
        case 4
            [Ks__ij,D1Ks__ijb,D2Ks__ijbg,D3Ks__ijbgd] = gaussianTKsC(q0.T,scales.^2,scaleweight.^2,int64(cdim),int64(L),int64(R));            
            D0Ks = tensor(Ks__ij,[L L],'ij');
            D1Ks = tensor(D1Ks__ijb,[L L cdim],'ijb');
            D2Ks = tensor(D2Ks__ijbg,[L L cdim cdim],'ijbg');
            D3Ks = tensor(D3Ks__ijbgd,[L L cdim cdim cdim],'ijbgd');
        otherwise
            assert(false);
    end
end

function M = dKs(x,y,dx,dy,r,sw)
    M = 2*fs.d1gamma(x,y,r,sw)*(x-y)'*(dx-dy);
end
function M = dN1Ks(x,y,dx,dy,r,sw)
    M = 4*fs.d2gamma(x,y,r,sw)*(x-y)'*(dx-dy)*(x-y)+2*fs.d1gamma(x,y,r,sw)*(dx-dy);
end
function M = dN2Ks(x,y,dx,dy,r,sw)
    M = -dN1Ks(x,y,dx,dy,r,sw);
end
function M = dD1N2Ks(x,y,dx,dy,r,sw)
    M = -4*fs.d2gamma(x,y,r,sw)*(x-y)'*(dx-dy)*eye(dim)...
        -8*fs.d3gamma(x,y,r,sw)*(x-y)'*(dx-dy)*(x-y)*(x-y)'...
        -4*fs.d2gamma(x,y,r,sw)*(dx-dy)*(x-y)'...
        -4*fs.d2gamma(x,y,r,sw)*(x-y)*(dx-dy)';
end
function M = dD2N2Ks(x,y,dx,dy,r,sw)
    M = -dD1N2Ks(x,y,dx,dy,r,sw);
end
function M = dD2D1N2Ksa(a,x,y,dx,dy,r,sw)
    M = 8*fs.d3gamma(x,y,r,sw)*(x-y)'*(dx-dy)*a*(x-y)'...
        +4*fs.d2gamma(x,y,r,sw)*a*(dx-dy)'...        
        +8*fs.d3gamma(x,y,r,sw)*(x-y)'*(dx-dy)*(x-y)*a'...        
        +4*fs.d2gamma(x,y,r,sw)*(dx-dy)*a'...        
        +8*fs.d3gamma(x,y,r,sw)*(x-y)'*(dx-dy)*(x-y)'*a*eye(dim)...        
        +4*fs.d2gamma(x,y,r,sw)*(dx-dy)'*a*eye(dim)...        
        +16*fs.d4gamma(x,y,r,sw)*(x-y)'*(dx-dy)*(x-y)'*a*(x-y)*(x-y)'...        
        +8*fs.d3gamma(x,y,r,sw)*(dx-dy)'*a*(x-y)*(x-y)'...        
        +8*fs.d3gamma(x,y,r,sw)*(x-y)'*a*(dx-dy)*(x-y)'...        
        +8*fs.d3gamma(x,y,r,sw)*(x-y)'*a*(x-y)*(dx-dy)';
end

fs.gamma = @(x,y,r,sw) 1/sw^2*exp(-sum((x-y).^2)/r^2);
fs.d1gamma = @(x,y,r,sw) -1/(sw^2*r^2)*exp(-sum((x-y).^2)/r^2);
fs.d2gamma = @(x,y,r,sw) 1/(sw^2*r^4)*exp(-sum((x-y).^2)/r^2);
fs.d3gamma = @(x,y,r,sw) -1/(sw^2*r^6)*exp(-sum((x-y).^2)/r^2);
fs.d4gamma = @(x,y,r,sw) 1/(sw^2*r^8)*exp(-sum((x-y).^2)/r^2);

fs.Ks = @Ks;
fs.N1Ks = @N1Ks;
fs.N2Ks = @N2Ks;
fs.D1N1Ks = @D1N1Ks;
fs.D2N1Ks = @D2N1Ks;
fs.D1N2Ks = @D1N2Ks;
fs.D2N2Ks = @D2N2Ks;
fs.D2D1N2Ksa = @D2D1N2Ksa;

fs.DaKs = @DaKs;

fs.TKs = @TKs;

fs.dKs = @dKs;
fs.dN1Ks = @dN1Ks;
fs.dN2Ks = @dN2Ks;
fs.dD1N2Ks = @dD1N2Ks;
fs.dD2N2Ks = @dD2N2Ks;
fs.dD2D1N2Ksa = @dD2D1N2Ksa;

end