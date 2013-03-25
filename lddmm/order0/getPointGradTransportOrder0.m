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

function gradTransport = getPointGradTransportOrder0(lddmmoptions)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim; % computations performed in cdim
L = lddmmoptions.L;
R = lddmmoptions.R;
CSP = lddmmoptions.CSP;
cCSP = lddmmoptions.cCSP;
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;
energyweight = lddmmoptions.energyweight;
epsilon = lddmmoptions.epsilon;

    function v0 = lgradTransport(v1, x, rhot)

        function dy = Gc(tt,ytt) % wrapper for cpu version of G
            t = intTime(tt,true,lddmmoptions);
            rhott = deval(rhot,t);
            if dim ~= cdim
                yt = rho2dTo3dOrder0(ytt,lddmmoptions);
            end    

            dy = fastPointGradTransportOrder0(yt,rhott,L,R,cdim,scales.^2,scaleweight.^2,energyweight);
            if dim ~= cdim
                dy = rho3dTo2dOrder0(dy,lddmmoptions);
            end
            
            dy = -intResult(dy,true,lddmmoptions); % sign for backwards integration already accounted for

            % debug
            if getOption(lddmmoptions,'testC')
                dy2 = G(tt,ytt);  
                assert(norm(dy-dy2) < epsilon);
            end
        end       
        [Ks D1Ks D2Ks] = gaussianKernels();
        function vt = Egradth(ttt)    
            vt = zeros(CSP*L,1);
            
            rhott = deval(rhot,ttt);
            if dim ~= cdim
                rhott = rho3dTo2dOrder0(rhott,lddmmoptions);
            end
            rhott = reshape(rhott,CSP,L);
            
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
        end
        function dy = G(tt,ytt)
            ytt = reshape(ytt,1,CSP*L);
            t = intTime(tt,true,lddmmoptions);          

            dy = zeros(size(ytt));

            rhott = deval(rhot,t);
            if dim ~= cdim
                rhott = rho3dTo2dOrder0(rhott,lddmmoptions);
            end
            rhott = reshape(rhott,CSP,L);

            for i = 1:L % particle
                xi = rhott(1:dim,i);
                dxi = ytt(1,CSP*(i-1)+(1:dim))';

                for l = 1:L % particle
                    xl = rhott(1:dim,l);   
                    dxl = ytt(1,CSP*(l-1)+(1:dim))';

                    ximxl = xi-xl;                    

                    for sl = 1:R % momentum
                        alsl = rhott(dim*(sl-1)+(1+dim:2*dim),l);
                        aisl = rhott(dim*(sl-1)+(1+dim:2*dim),i);

                        ksl = Ks(xi,xl,scales(sl),scaleweight(sl));
                        d1ksl = D1Ks(xi,xl,scales(sl),scaleweight(sl));
                        d1kslXximxl = d1ksl*ximxl;
                        d2ksl = D2Ks(xi,xl,scales(sl),scaleweight(sl));                  

                        dalsl = ytt(1,CSP*(l-1)+dim*(sl-1)+(1+dim:2*dim))';

                        % dx
                        dy(1,CSP*(i-1)+(1:dim)) = dy(1,CSP*(i-1)+(1:dim)) ...
                            + (2*d1kslXximxl*(aisl'*dxl+alsl'*dxi))'; 

                        for si = 1:R % scale                
                            aisi = rhott(dim*(si-1)+(1+dim:2*dim),i);
                            alsi = rhott(dim*(si-1)+(1+dim:2*dim),l);
                            daisi = ytt(1,CSP*(i-1)+dim*(si-1)+(1+dim:2*dim))';
                            daisl = ytt(1,CSP*(i-1)+dim*(sl-1)+(1+dim:2*dim))';

                            ksi = Ks(xi,xl,scales(si),scaleweight(si));
                            d1ksi = D1Ks(xi,xl,scales(si),scaleweight(si));
                            d2ksi = D2Ks(xi,xl,scales(si),scaleweight(si));                               

                            % dx
                            dy(1,CSP*(i-1)+(1:dim)) = dy(1,CSP*(i-1)+(1:dim)) ...
                                + (-2*d1ksi*(aisl'*alsi*daisl-aisi'*alsl*dalsl) ...
                                - 4*d2ksi*ximxl'*(aisl'*alsi*daisl-aisi'*alsl*dalsl)*ximxl)';  

                            % da
                            dy(1,CSP*(i-1)+dim*(si-1)+(1+dim:2*dim)) = dy(1,CSP*(i-1)+dim*(si-1)+(1+dim:2*dim)) ...
                                + (1/R*ksi*dxl ...
                                - 2*ximxl'*(d1ksl*daisi-d1ksi*dalsl)*alsl)';
                        end
                    end
                end
            end

            dy = dy + energyweight(1)*Egradth(t)';
            dy = -intResult(dy,true,lddmmoptions); % sign for backwards integration already accounted for            
            
            dy = reshape(dy,CSP*L,1);
        end

        % integrate
        options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        vt = ode45(@Gc,[0 1],reshape(v1,CSP,L),options); % solve backwards, cpu
%         vt = ode45(@G,[0 1],v1,options); % solve backwards, slooow matlab
        assert(vt.x(end) == 1);
        v0 = reshape(deval(vt,1),CSP*L,1);
        
    end

gradTransport = @lgradTransport;

end