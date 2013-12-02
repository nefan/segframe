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

function pathEnergy = getPointPathEnergyOrder0(lddmmoptions)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim; % computations performed in cdim
L = lddmmoptions.L;
R = lddmmoptions.R;
CSP = lddmmoptions.CSP;
cCSP = lddmmoptions.cCSP;
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;

[Ks D1Ks D2Ks] = gaussianKernels();

    function [E v] = lpathEnergy(x,rhot,varargin)
        tend = 1;
        if size(varargin,2) > 0
            tend = varargin{1};
        end
        
        function vt = gradEt(tt,y)
            t = intTime(tt,true,lddmmoptions);            
            vt = zeros(cCSP*L,1);
            rhott = reshape(deval(rhot,t),cCSP,L);
            for i = 1:L % particle
                xi = rhott(1:dim,i);

                for ll = 1:L % particle
                    xl = rhott(1:dim,ll);

                    for ss = 1:R                                           
                        al = rhott(dim*(ss-1)+(1+dim:2*dim),ll); 
                        ai = rhott(dim*(ss-1)+(1+dim:2*dim),i);                                                                                                               

                        ks = Ks(xi,xl,scales(ss),scaleweight(ss));
                        d1ks = D1Ks(xi,xl,scales(ss),scaleweight(ss));

                        vt(CSP*(i-1)+dim*(ss-1)+(1+dim:2*dim)) = vt(CSP*(i-1)+dim*(ss-1)+(1+dim:2*dim)) + 2*ks*al;
                        vt(CSP*(i-1)+(1:dim)) = vt(CSP*(i-1)+(1:dim)) + 4*d1ks*ai'*al*(xi-xl);
                    end
                end
            end
            
            vt = -intResult(vt,true,lddmmoptions); % sign for backwards integration already accounted for
        end
        
        function Et = Gc(tt) % wrapper for C version of G
            t = intTime(tt,false,lddmmoptions);
            rhott = deval(rhot,t);         

            Et = fastPointPathEnergyOrder0(rhott,L,R,cdim,scales.^2,scaleweight.^2);

            % debug
%             [t Et]
            if getOption(lddmmoptions,'testC')
                Et2 = G(tt);  
                assert(norm(Et2-Et) < 10e-12);
            end
        end

        function Et = G(tt)
            t = intTime(tt,false,lddmmoptions);
            rhott = deval(rhot,t);
            if dim ~= cdim
                rhott = rho3dTo2dOrder0(rhott,lddmmoptions);
            end
            rhott = reshape(rhott,CSP,L);

            Et = 0;

            for s = 1:R
                for i = 1:L
                    for l = 1:L                
                        Et = Et + rhott(dim*(s-1)+(1+dim:2*dim),i)'*Ks(rhott(1:dim,i),rhott(1:dim,l),scales(s),scaleweight(s))*rhott(dim*(s-1)+(1+dim:2*dim),l);
                    end
                end
            end
        end

        E = integrate(@Gc,0,tend); % fast C version
        % E = energyweight(1)*integrate(@G,0,tend); % sloow matlab version
        assert(E >= 0);
        
        if nargout > 1
            v1 = zeros(CSP*L,1);
            options = odeset('RelTol',1e-6,'AbsTol',1e-6);
            vt = ode45(@gradEt,[0 tend],v1,options);        
            v = deval(vt,0);
        end
    end

pathEnergy = @lpathEnergy;

end