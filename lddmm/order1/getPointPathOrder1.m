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
iCSP = lddmmoptions.iCSP;
fCSP = lddmmoptions.fCSP;
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;
epsilon = lddmmoptions.epsilon;

% indexing functions
iImu = lddmmoptions.iImu;
iImuj = lddmmoptions.iImuj;
fIphi = lddmmoptions.fIphi;
fIDphi = lddmmoptions.fIDphi;
fIDphic = lddmmoptions.fIDphic;
fImu = lddmmoptions.fImu;

ks = dkernelsGaussian(cdim);

    function Gt = pointPathOrder1(x)
        x = reshape(x,iCSP,L);
        rhoj0 = zeros(R*dim^2,L);
        for s=1:R
            rhoj0(dim^2*(s-1)+(1:dim^2),:) = x(iImuj(s),:);
        end
        if dim ~= cdim
            rhoj0 = reshape(rhoj2dTo3dOrder1(rhoj0,lddmmoptions),R*cdim^2,L);
        end        
        
        function dGt = Gc(t,Gt) % wrapper for cpu version of G
            assert(false);

            dGt = fastPointPathOrder1(t,Gt,reshape(rhoj0,[],1),L,R,cdim,scales.^2,scaleweight.^2);

            % debug
            if getOption(lddmmoptions,'testC')
                dGt2 = G(t,Gt);  
%                 norm(dGt-dGt2)
                assert(norm(dGt-dGt2) < epsilon);
            end
        end           
    
        function dGt = G(t,Gt)  % slooow version
            Gt = reshape(Gt,fCSP,L);

            dGt = zeros(size(Gt));

            for n = 1:L % particle                
                xn = Gt(fIphi(),n);
                Dphin = reshape(Gt(fIDphi(),n),cdim,cdim);

                for i = 1:L % particle                    
                    xi = Gt(1:cdim,i);
                    Dphii = reshape(Gt(fIDphi(),i),cdim,cdim);            

                    for si = 1:R % scale
                        mui = Gt(fImu(si),i);

                        % position
                        dGt(fIphi(),n) = dGt(fIphi(),n)+...
                            ks.Ks(xi,xn,scales(si),scaleweight(si))*mui;
                        for j = 1:cdim
                            muij = Dphii'\rhoj0(cdim^2*(si-1)+cdim*(j-1)+(1:cdim),i);
                            dGt(fIphi(),n) = dGt(fIphi(),n)...
                                +ks.N1Ks(xi,xn,scales(si),scaleweight(si))'*Dphii(:,j)*muij;
                        end
                        % dPhi
                        for c = 1:cdim
                            dGt(fIDphic(c),n) = dGt(fIDphic(c),n)...
                                +ks.N2Ks(xi,xn,scales(si),scaleweight(si))'*Dphin(:,c)*mui;
                            for j = 1:cdim
                                muij = Dphii'\rhoj0(cdim^2*(si-1)+cdim*(j-1)+(1:cdim),i);
                                dGt(fIDphic(c),n) = dGt(fIDphic(c),n)...
                                    +(ks.D1N2Ks(xi,xn,scales(si),scaleweight(si))*Dphii(:,j))'*Dphin(:,c)*muij;
                            end                
                        end

                        for sn = 1:R % scale
                            mun = Gt(fImu(sn),n);

                            % mu
                            dGt(fImu(sn),n) = dGt(fImu(sn),n)...
                                -mun'*mui*ks.N2Ks(xi,xn,scales(si),scaleweight(si));
                            for j = 1:cdim
                                muij = Dphii'\rhoj0(cdim^2*(si-1)+cdim*(j-1)+(1:cdim),i);
                                munj = Dphin'\rhoj0(cdim^2*(sn-1)+cdim*(j-1)+(1:cdim),n);
                                dGt(fImu(sn),n) = dGt(fImu(sn),n)...
                                    -munj'*mui*ks.D2N2Ks(xi,xn,scales(si),scaleweight(si))*Dphin(:,j) ...
                                    -mun'*muij*ks.D1N2Ks(xi,xn,scales(si),scaleweight(si))*Dphii(:,j);
                                    %-(munj'*mui-mun'*muij)*D2N2Ks(xi,xn,scales(si),scaleweight(si))*Dphin(:,j);                       
                                for jm = 1:cdim
                                    munjm = Dphin'\rhoj0(cdim^2*(sn-1)+cdim*(jm-1)+(1:cdim),n);
                                    dGt(fImu(sn),n) = dGt(fImu(sn),n)...
                                        -munjm'*muij*ks.D2D1N2Ksa(Dphii(:,j),xi,xn,scales(si),scaleweight(si))*Dphin(:,jm);
                                end                                                            
                            end                                    
                        end                
                    end
                end
            end

            dGt = reshape(dGt,L*fCSP,1);
        end

        options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        initial = zeros(fCSP,L);
        initial(fIphi(),:) = [moving; zeros(cdim-dim,L)]; % position
        initial(fIDphi(),:) = repmat(reshape(eye(cdim),cdim^2,1),1,L); % Dphi
        for s=1:R
            initial(fImu(s),:) = [x(iImu(s),:); zeros(cdim-dim,L)]; % mu/rho0
        end
%         Gt = ode45(@Gc,[0 1],reshape(initial,L*fCSP,1),options); % fast native version
        Gt = ode45(@G,[0 1],reshape(initial,L*fCSP,1),options); % matlab version
        assert(Gt.x(end) == 1); % if not, integration failed

    end

pointPath = @pointPathOrder1;

end