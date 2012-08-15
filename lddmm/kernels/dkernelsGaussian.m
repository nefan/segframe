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
fs.D2N1Ks = @D2N1Ks;
fs.D1N2Ks = @D1N2Ks;
fs.D2N2Ks = @D2N2Ks;
fs.D2D1N2Ksa = @D2D1N2Ksa;

fs.dKs = @dKs;
fs.dN1Ks = @dN1Ks;
fs.dN2Ks = @dN2Ks;
fs.dD1N2Ks = @dD1N2Ks;
fs.dD2N2Ks = @dD2N2Ks;
fs.dD2D1N2Ksa = @dD2D1N2Ksa;

end