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

function gradTransport = getPointGradTransportOrder1(lddmmoptions)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim;
L = lddmmoptions.L;
R = lddmmoptions.R;
CSP = lddmmoptions.CSP;
cCSP = lddmmoptions.cCSP;
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;
energyweight = lddmmoptions.energyweight;
epsilon = lddmmoptions.epsilon;

phioffset = 0;
Dphioffset = phioffset+cdim;
muoffset = Dphioffset+cdim*cdim;
mujoffset = muoffset+cdim;
gradCSP = 2*dim+2*dim^2;
cgradCSP = 2*cdim+2*cdim^2;

ks = dkernelsGaussian(dim);

    function v0 = lgradTransport(v1, x, Gt)
        
        x = reshape(x,(1+dim)*dim,L);
        rho0 = x(1:dim,:);
        rhoj0 = x(dim+(1:dim^2),:);   
        if dim ~= cdim
            rho0 = reshape(rho2dTo3dOrder1(rho0,lddmmoptions),cdim,L);
            rhoj0 = reshape(rhoj2dTo3dOrder1(rhoj0,lddmmoptions),cdim^2,L);
        end
        M = dkernelDerivativeMatrix(rho0,rhoj0,lddmmoptions);

%         function v = Egrad()    
%             v = zeros(cdim*(1+cdim),L);
%             G0 = reshape(deval(Gt,0),cdim*(1+cdim+1),L);
%             for n = 1:L
%                 xn = G0(1:cdim,n);
%                 an = G0(cdim*(1+cdim)+(1:cdim),n); 
% 
%                 for ll = 1:L % particle
%                     xl = G0(1:cdim,ll);
%                     al = G0(cdim*(1+cdim)+(1:cdim),ll);  
%                     Ks = ks.Ks(xl,xn,scales(s),scaleweight(s));
%                     N1 = ks.N1Ks(xl,xn,scales(s),scaleweight(s));
%                     D2N1 = ks.D2N1Ks(xl,xn,scales(s),scaleweight(s));
% 
%                     % an
%                     v(1:cdim,n) = v(1:cdim,n) + 2*Ks*al;
% 
%                     for j = 1:cdim
%                         alj = rhoj0(cdim*(j-1)+(1:cdim),ll);             
%                         % an
%                         v(1:cdim,n) = v(1:cdim,n) + 2*N1(j)*alj;
% 
%                         % anj
%                         v(cdim+cdim*(j-1)+(1:cdim),n) = v(cdim+cdim*(j-1)+(1:cdim),n) - 2*N1(j)*an;
% 
%                         for jm = 1:cdim
%                             aljm = rhoj0(cdim*(j-1)+(1:cdim),ll);             
% 
%                             % anj
%                             v(cdim+cdim*(j-1)+(1:cdim),n) = v(cdim+cdim*(j-1)+(1:cdim),n) + 2*D2N1(jm,j)*aljm;
%                         end
%                     end            
%                 end
%             end
% 
%             v = reshape(v,cdim*(1+cdim)*L,1);
%         end
        
        function dy = Gc(tt,ytt) % wrapper for cpu version of G
            t = intTime(tt,true,lddmmoptions);
            Gtt = deval(Gt,t);
            
            dy = fastPointGradTransportOrder1(ytt,Gtt,reshape(rhoj0,[],1),L,R,cdim,scales.^2,scaleweight.^2,energyweight);
%             dy = dy + Gtest(t,ytt); % until everything is implemented

            dy = -intResult(dy,true,lddmmoptions); % sign for backwards integration already accounted for

            % debug
            if getOption(lddmmoptions,'testC')
                dy2 = G(tt,ytt);  
%                 norm(dy-dy2)
                assert(norm(dy-dy2) < epsilon);
            end
        end   
        
        function dy = Gtest(tt,ytt)
            ytt = reshape(ytt,cgradCSP,L);
            t = intTime(t,false,lddmmoptions);

            dy = zeros(size(ytt));

            Gtt = reshape(deval(Gt,t),dim*(1+dim+1),L);
            if dim ~= cdim
                Gtt = reshape(Gt2dTo3dOrder1(Gtt,lddmmoptions),cCSP,L);
            end 

            for n = 1:L % update all particles
                for i = 1:L % multiply
%                     dy(phioffset+(1:cdim),n) = dy(phioffset+(1:cdim),n) ...
%                         + M.aphiphi(i,n,Gtt)'*ytt(phioffset+(1:cdim),i) ...
%                         + M.amuphi(i,n,Gtt)'*ytt(muoffset+(1:cdim),i);
%                     for ll = 1:cdim
%                         dy(phioffset+(1:cdim),n) = dy(phioffset+(1:cdim),n) ...
%                             + M.aDphiphi(i,n,ll,Gtt)'*ytt(Dphioffset+cdim*(ll-1)+(1:cdim),i);
%                     end            
%                     for jj = 1:cdim
%                         dy(phioffset+(1:cdim),n) = dy(phioffset+(1:cdim),n) ...
%                             + M.amujphi(i,n,jj,Gtt)'*ytt(mujoffset+cdim*(jj-1)+(1:cdim),i);                
%                     end

%                     for kk = 1:cdim
%                         dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) = dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) ...
%                             + M.aphiDphi(i,n,kk,Gtt)'*ytt(phioffset+(1:cdim),i) ...
%                             + M.amuDphi(i,n,kk,Gtt)'*ytt(muoffset+(1:cdim),i);
%                         for ll = 1:cdim
%                             dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) = dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) ...
%                                 + M.aDphiDphi(i,n,ll,kk,Gtt)'*ytt(Dphioffset+cdim*(ll-1)+(1:cdim),i);
%                         end            
%                         for jj = 1:cdim
%                             dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) = dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) ...
%                                 + M.amujDphi(i,n,jj,kk,Gtt)'*ytt(mujoffset+cdim*(jj-1)+(1:cdim),i);                
%                         end           
%                     end

%                     dy(muoffset+(1:cdim),n) = dy(muoffset+(1:cdim),n) ...
%                         + M.aphimu(i,n,Gtt)'*ytt(phioffset+(1:cdim),i) ...
%                         + M.amumu(i,n,Gtt)'*ytt(muoffset+(1:cdim),i);
%                     for ll = 1:cdim
%                         dy(muoffset+(1:cdim),n) = dy(muoffset+(1:cdim),n) ...
%                             + M.aDphimu(i,n,ll,Gtt)'*ytt(Dphioffset+cdim*(ll-1)+(1:cdim),i);
%                     end            
%                     for jj = 1:cdim
%                         dy(muoffset+(1:cdim),n) = dy(muoffset+(1:cdim),n) ...
%                             + M.amujmu(i,n,jj,Gtt)'*ytt(mujoffset+cdim*(jj-1)+(1:cdim),i);                
%                     end

%                     for jjm = 1:cdim
%                         dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) = dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) ...
%                             + M.aphimuj(i,n,jjm,Gtt)'*ytt(phioffset+(1:cdim),i) ...
%                             + M.amumuj(i,n,jjm,Gtt)'*ytt(muoffset+(1:cdim),i);
%                         for ll = 1:cdim
%                             dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) = dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) ...
%                                 + M.aDphimuj(i,n,ll,jjm,Gtt)'*ytt(Dphioffset+cdim*(ll-1)+(1:cdim),i);
%                         end            
%                         for jj = 1:cdim
%                             dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) = dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) ...
%                                 + M.amujmuj(i,n,jj,jjm,Gtt)'*ytt(mujoffset+cdim*(jj-1)+(1:cdim),i);                
%                         end           
%                     end            
                end
            end

            dy = reshape(dy,cgradCSP*L,1);
        end        

        function dy = G(tt,ytt)
            ytt = reshape(ytt,cgradCSP,L);
            t = intTime(tt,true,lddmmoptions);

            dy = zeros(size(ytt));

            Gtt = reshape(deval(Gt,t),cdim*(1+cdim+1),L);

            for n = 1:L % update all particles
                for i = 1:L % multiply
                    dy(phioffset+(1:cdim),n) = dy(phioffset+(1:cdim),n) ...
                        + M.aphiphi(i,n,Gtt)'*ytt(phioffset+(1:cdim),i) ...
                        + M.amuphi(i,n,Gtt)'*ytt(muoffset+(1:cdim),i);
                    for ll = 1:cdim
                        dy(phioffset+(1:cdim),n) = dy(phioffset+(1:cdim),n) ...
                            + M.aDphiphi(i,n,ll,Gtt)'*ytt(Dphioffset+cdim*(ll-1)+(1:cdim),i);
                    end            
                    for jj = 1:cdim
                        dy(phioffset+(1:cdim),n) = dy(phioffset+(1:cdim),n) ...
                            + M.amujphi(i,n,jj,Gtt)'*ytt(mujoffset+cdim*(jj-1)+(1:cdim),i);                
                    end

                    for kk = 1:cdim
                        dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) = dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) ...
                            + M.aphiDphi(i,n,kk,Gtt)'*ytt(phioffset+(1:cdim),i) ...
                            + M.amuDphi(i,n,kk,Gtt)'*ytt(muoffset+(1:cdim),i);
                        for ll = 1:cdim
                            dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) = dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) ...
                                + M.aDphiDphi(i,n,ll,kk,Gtt)'*ytt(Dphioffset+cdim*(ll-1)+(1:cdim),i);
                        end            
                        for jj = 1:cdim
                            dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) = dy(Dphioffset+cdim*(kk-1)+(1:cdim),n) ...
                                + M.amujDphi(i,n,jj,kk,Gtt)'*ytt(mujoffset+cdim*(jj-1)+(1:cdim),i);                
                        end           
                    end

                    dy(muoffset+(1:cdim),n) = dy(muoffset+(1:cdim),n) ...
                        + M.aphimu(i,n,Gtt)'*ytt(phioffset+(1:cdim),i) ...
                        + M.amumu(i,n,Gtt)'*ytt(muoffset+(1:cdim),i);
                    for ll = 1:cdim
                        dy(muoffset+(1:cdim),n) = dy(muoffset+(1:cdim),n) ...
                            + M.aDphimu(i,n,ll,Gtt)'*ytt(Dphioffset+cdim*(ll-1)+(1:cdim),i);
                    end            
                    for jj = 1:cdim
                        dy(muoffset+(1:cdim),n) = dy(muoffset+(1:cdim),n) ...
                            + M.amujmu(i,n,jj,Gtt)'*ytt(mujoffset+cdim*(jj-1)+(1:cdim),i);                
                    end

                    for jjm = 1:cdim
                        dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) = dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) ...
                            + M.aphimuj(i,n,jjm,Gtt)'*ytt(phioffset+(1:cdim),i) ...
                            + M.amumuj(i,n,jjm,Gtt)'*ytt(muoffset+(1:cdim),i);
                        for ll = 1:cdim
                            dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) = dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) ...
                                + M.aDphimuj(i,n,ll,jjm,Gtt)'*ytt(Dphioffset+cdim*(ll-1)+(1:cdim),i);
                        end            
                        for jj = 1:cdim
                            dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) = dy(mujoffset+cdim*(jjm-1)+(1:cdim),n) ...
                                + M.amujmuj(i,n,jj,jjm,Gtt)'*ytt(mujoffset+cdim*(jj-1)+(1:cdim),i);                
                        end           
                    end            
                end
            end

            dy = -intResult(dy,true,lddmmoptions); % sign for backwards integration already accounted for
            dy = reshape(dy,cgradCSP*L,1);
        end

        % function dy = Gforward(t,ytt)
        %     dy1 = Gforward1(t,ytt);
        %     dy2 = Gforward2(t,ytt);
        %     
        %     e = norm(dy1-dy2);
        %     if e > 1e-5
        %         reshape(dy1-dy2,[],L)
        %         e
        %     end
        % %     assert(e < 1e-5);
        %     
        %     dy = dy1;
        % end
        % function dy = Gforward1(t,ytt)
        %     ytt = reshape(ytt,gradCSP,L);
        %     
        %     dy = zeros(size(ytt));
        %     
        %     Gtt = reshape(deval(Gt,t),dim*(1+dim+1),L);
        %     
        %     for n = 1:L % update all particles
        %         for i = 1:L % multiply
        %             dy(phioffset+(1:dim),n) = dy(phioffset+(1:dim),n) ...
        %                 + M.aphiphi(n,i,Gtt)*ytt(phioffset+(1:dim),i) ...
        %                 + M.aphimu(n,i,Gtt)*ytt(muoffset+(1:dim),i);
        %             for kk = 1:dim
        %                 dy(phioffset+(1:dim),n) = dy(phioffset+(1:dim),n) ...
        %                     + M.aphiDphi(n,i,kk,Gtt)*ytt(Dphioffset+dim*(kk-1)+(1:dim),i);
        %             end            
        %             for jj = 1:dim
        %                 dy(phioffset+(1:dim),n) = dy(phioffset+(1:dim),n) ...
        %                     + M.aphimuj(n,i,jj,Gtt)*ytt(mujoffset+dim*(jj-1)+(1:dim),i);                
        %             end
        % 
        %             for ll = 1:dim
        %                 dy(Dphioffset+dim*(ll-1)+(1:dim),n) = dy(Dphioffset+dim*(ll-1)+(1:dim),n) ...
        %                     + M.aDphiphi(n,i,ll,Gtt)*ytt(phioffset+(1:dim),i) ...
        %                     + M.aDphimu(n,i,ll,Gtt)*ytt(muoffset+(1:dim),i);
        %                 for kk = 1:dim
        %                     dy(Dphioffset+dim*(ll-1)+(1:dim),n) = dy(Dphioffset+dim*(ll-1)+(1:dim),n) ...
        %                         + M.aDphiDphi(n,i,ll,kk,Gtt)*ytt(Dphioffset+dim*(kk-1)+(1:dim),i);
        %                 end            
        %                 for jj = 1:dim
        %                     dy(Dphioffset+dim*(ll-1)+(1:dim),n) = dy(Dphioffset+dim*(ll-1)+(1:dim),n) ...
        %                         + M.aDphimuj(n,i,ll,jj,Gtt)*ytt(mujoffset+dim*(jj-1)+(1:dim),i);                
        %                 end           
        %             end
        %             
        %             dy(muoffset+(1:dim),n) = dy(muoffset+(1:dim),n) ...
        %                 + M.amuphi(n,i,Gtt)*ytt(phioffset+(1:dim),i) ...
        %                 + M.amumu(n,i,Gtt)*ytt(muoffset+(1:dim),i);
        %             for kk = 1:dim
        %                 dy(muoffset+(1:dim),n) = dy(muoffset+(1:dim),n) ...
        %                     + M.amuDphi(n,i,kk,Gtt)*ytt(Dphioffset+dim*(kk-1)+(1:dim),i);
        %             end            
        %             for jj = 1:dim
        %                 dy(muoffset+(1:dim),n) = dy(muoffset+(1:dim),n) ...
        %                     + M.amumuj(n,i,jj,Gtt)*ytt(mujoffset+dim*(jj-1)+(1:dim),i);                
        %             end
        % 
        %             for jj = 1:dim
        %                 dy(mujoffset+dim*(jj-1)+(1:dim),n) = dy(mujoffset+dim*(jj-1)+(1:dim),n) ...
        %                     + M.amujphi(n,i,jj,Gtt)*ytt(phioffset+(1:dim),i) ...
        %                     + M.amujmu(n,i,jj,Gtt)*ytt(muoffset+(1:dim),i);
        %                 for kk = 1:dim
        %                     dy(mujoffset+dim*(jj-1)+(1:dim),n) = dy(mujoffset+dim*(jj-1)+(1:dim),n) ...
        %                         + M.amujDphi(n,i,jj,kk,Gtt)*ytt(Dphioffset+dim*(kk-1)+(1:dim),i);
        %                 end            
        %                 for jjm = 1:dim
        %                     dy(mujoffset+dim*(jj-1)+(1:dim),n) = dy(mujoffset+dim*(jj-1)+(1:dim),n) ...
        %                         + M.amujmuj(n,i,jj,jjm,Gtt)*ytt(mujoffset+dim*(jjm-1)+(1:dim),i);                
        %                 end           
        %             end            
        %         end
        %     end
        % 
        %     dy = reshape(dy,gradCSP*L,1);
        % end
        % function dy = Gforward2(t,ytt)
        %     ytt = reshape(ytt,gradCSP,L);
        %     
        %     dy = zeros(size(ytt));
        %     
        %     Gtt = reshape(deval(Gt,t),dim*(1+dim+1),L);
        %     
        %     for n = 1:L % update all particles
        %         xn = M.x(n,Gtt);
        %         dxn = ytt(phioffset+(1:dim),n);
        %         dmun = ytt(muoffset+(1:dim),n);
        %         
        %         for i = 1:L % multiply
        %             xi = M.x(i,Gtt);
        %             dxi = ytt(phioffset+(1:dim),i);
        %             dmui = ytt(muoffset+(1:dim),i);
        %             
        %             dy(phioffset+(1:dim),n) = dy(phioffset+(1:dim),n) ...
        %                 + ks.dKs(xi,xn,dxi,dxn,scales(s),scaleweight(s))*M.mu(i,Gtt)+ks.gamma(xi,xn,scales(s),scaleweight(s))*dmui;
        %             
        %             for jj = 1:dim
        %                 dDphiij = ytt(Dphioffset+dim*(jj-1)+(1:dim),i);
        %                 dmuij = ytt(mujoffset+dim*(jj-1)+(1:dim),i);
        %                 
        %                 dy(phioffset+(1:dim),n) = dy(phioffset+(1:dim),n) ...
        %                     + ks.dN1Ks(xi,xn,dxi,dxn,scales(s),scaleweight(s))'*M.Dphil(i,jj,Gtt)*M.muj(i,jj,Gtt)...
        %                     + ks.N1Ks(xi,xn,scales(s),scaleweight(s))'*dDphiij*M.muj(i,jj,Gtt)...    
        %                     + ks.N1Ks(xi,xn,scales(s),scaleweight(s))'*M.Dphil(i,jj,Gtt)*dmuij;                
        %             end
        % 
        %             for ll = 1:dim
        %                 dDphinl = ytt(Dphioffset+dim*(ll-1)+(1:dim),n);
        %                 
        %                 dy(Dphioffset+dim*(ll-1)+(1:dim),n) = dy(Dphioffset+dim*(ll-1)+(1:dim),n) ...
        %                     + ks.dN2Ks(xi,xn,dxi,dxn,scales(s),scaleweight(s))'*M.Dphil(n,ll,Gtt)*M.mu(i,Gtt)...
        %                     + ks.N2Ks(xi,xn,scales(s),scaleweight(s))'*dDphinl*M.mu(i,Gtt)...
        %                     + ks.N2Ks(xi,xn,scales(s),scaleweight(s))'*M.Dphil(n,ll,Gtt)*dmui;
        %                 
        %                 for jj = 1:dim
        %                     dDphiij = ytt(Dphioffset+dim*(jj-1)+(1:dim),i);
        %                     dmuij = ytt(mujoffset+dim*(jj-1)+(1:dim),i);                    
        %                     
        %                     dy(Dphioffset+dim*(ll-1)+(1:dim),n) = dy(Dphioffset+dim*(ll-1)+(1:dim),n) ...
        %                         + (ks.dD1N2Ks(xi,xn,dxi,dxn,scales(s),scaleweight(s))*M.Dphil(i,jj,Gtt))'*M.Dphil(n,ll,Gtt)*M.muj(i,jj,Gtt)...
        %                         + (ks.D1N2Ks(xi,xn,scales(s),scaleweight(s))*dDphiij)'*M.Dphil(n,ll,Gtt)*M.muj(i,jj,Gtt)...
        %                         + (ks.D1N2Ks(xi,xn,scales(s),scaleweight(s))*M.Dphil(i,jj,Gtt))'*dDphinl*M.muj(i,jj,Gtt)...
        %                         + (ks.D1N2Ks(xi,xn,scales(s),scaleweight(s))*M.Dphil(i,jj,Gtt))'*M.Dphil(n,ll,Gtt)*dmuij;
        %                 end           
        %             end
        %             
        %             dy(muoffset+(1:dim),n) = dy(muoffset+(1:dim),n) ...
        %                 - (dmun'*M.mu(i,Gtt)+M.mu(n,Gtt)'*dmui)*ks.N2Ks(xi,xn,scales(s),scaleweight(s))...
        %                 - M.mu(n,Gtt)'*M.mu(i,Gtt)*ks.dN2Ks(xi,xn,dxi,dxn,scales(s),scaleweight(s));
        %             for jj = 1:dim
        %                 dmunj = ytt(mujoffset+dim*(jj-1)+(1:dim),n);
        %                 dmuij = ytt(mujoffset+dim*(jj-1)+(1:dim),i);
        %                 dDphiij = ytt(Dphioffset+dim*(jj-1)+(1:dim),i);
        %                 dDphinj = ytt(Dphioffset+dim*(jj-1)+(1:dim),n);
        %                 
        %                 dy(muoffset+(1:dim),n) = dy(muoffset+(1:dim),n) ...
        %                     - (dmunj'*M.mu(i,Gtt)+M.muj(n,jj,Gtt)'*dmui-dmun'*M.muj(i,jj,Gtt)-M.mu(n,Gtt)'*dmuij)...
        %                        *ks.D2N2Ks(xi,xn,scales(s),scaleweight(s))*M.Dphil(n,jj,Gtt)...
        %                     - (M.muj(n,jj,Gtt)'*M.mu(i,Gtt)-M.mu(n,Gtt)'*M.muj(i,jj,Gtt))...
        %                        *ks.dD2N2Ks(xi,xn,dxi,dxn,scales(s),scaleweight(s))*M.Dphil(n,jj,Gtt)...
        %                     - (M.muj(n,jj,Gtt)'*M.mu(i,Gtt)-M.mu(n,Gtt)'*M.muj(i,jj,Gtt))...
        %                        *ks.D2N2Ks(xi,xn,scales(s),scaleweight(s))*dDphinj;  
        %                    
        %                    for jjm = 1:dim
        %                        dmunjm = ytt(mujoffset+dim*(jjm-1)+(1:dim),n);
        %                        dDphinjm = ytt(Dphioffset+dim*(jjm-1)+(1:dim),n);
        %                        
        %                        dy(muoffset+(1:dim),n) = dy(muoffset+(1:dim),n) ...
        %                             - (dmunjm'*M.muj(i,jj,Gtt)+M.muj(n,jjm,Gtt)'*dmuij)...
        %                                *ks.D2D1N2Ksa(M.Dphil(i,jj,Gtt),xi,xn,scales(s),scaleweight(s))*M.Dphil(n,jjm,Gtt)...     
        %                             - (M.muj(n,jjm,Gtt)'*M.muj(i,jj,Gtt))...
        %                                *ks.D2D1N2Ksa(dDphiij,xi,xn,scales(s),scaleweight(s))*M.Dphil(n,jjm,Gtt)...
        %                             - (M.muj(n,jjm,Gtt)'*M.muj(i,jj,Gtt))...
        %                                *ks.dD2D1N2Ksa(M.Dphil(i,jj,Gtt),xi,xn,dxi,dxn,scales(s),scaleweight(s))*M.Dphil(n,jjm,Gtt)...                               
        %                             - (M.muj(n,jjm,Gtt)'*M.muj(i,jj,Gtt))...
        %                                *ks.D2D1N2Ksa(M.Dphil(i,jj,Gtt),xi,xn,scales(s),scaleweight(s))*dDphinjm;                               
        %                    end
        %             end
        %         end
        %     end
        % 
        %     for n = 1:L % update all particles
        %         xn = M.x(n,Gtt);
        %         dxn = ytt(phioffset+(1:dim),n);
        %         dmun = ytt(muoffset+(1:dim),n);
        %         
        %         for i = 1:L % multiply
        %             xi = M.x(i,Gtt);
        %             dxi = ytt(phioffset+(1:dim),i);
        %             dmui = ytt(muoffset+(1:dim),i);   
        %             
        %             if i == n
        %                 for jj = 1:dim
        %                     
        %                     % debug
        %                     dDphin = reshape(ytt(Dphioffset+(1:dim*dim),n),dim,dim);
        %                     a0nj = rhoj0(dim*(jj-1)+(1:dim),n);
        %                     da0nj = da0j(dim*(jj-1)+(1:dim),n);
        % %                     dmunjCheck = (inv(M.Dphi(n,Gtt))*dDphin*inv(M.Dphi(n,Gtt)))'*a0nj+inv(M.Dphi(n,Gtt))'*da0nj;
        %                     dmunjCheck = -M.Dphi(n,Gtt)'\(dDphin'*(M.Dphi(n,Gtt)'\a0nj))+M.Dphi(n,Gtt)'\da0nj;
        %                     dmunj = ytt(mujoffset+dim*(jj-1)+(1:dim),n);
        %                     e = norm(dmunj-dmunjCheck);
        %                     if e > 1e-4
        % %                         norm(inv(M.Dphi(n,Gtt))*M.Dphi(n,Gtt)-eye(dim))
        % %                         norm(M.Dphi(n,Gtt)*inv(M.Dphi(n,Gtt))-eye(dim))
        %                         fprintf('dmunj error (t=%f): %f\n',t,e);
        %                     end
        % %                     assert(e < 1e-6);
        % 
        %                     dy(mujoffset+dim*(jj-1)+(1:dim),n) = dy(mujoffset+dim*(jj-1)+(1:dim),n) ...
        %                             -M.Dphi(n,Gtt)'\(M.dtDphi(n,Gtt)'*dmunj);
        %                     for jjm = 1:dim
        %                         dDphinjm = ytt(Dphioffset+dim*(jjm-1)+(1:dim),n);
        %                         dtdDphinjm = dy(Dphioffset+dim*(jjm-1)+(1:dim),n);
        %                         ejm = zeros(dim,1); ejm(jjm) = 1;
        % 
        %                         dy(mujoffset+dim*(jj-1)+(1:dim),n) = dy(mujoffset+dim*(jj-1)+(1:dim),n) ...
        %                             +(M.Dphi(n,Gtt)'\ejm)*(inv(M.Dphi(n,Gtt))'*M.dtDphi(n,Gtt)'*M.muj(n,jj,Gtt))'*dDphinjm...
        %                             -(M.Dphi(n,Gtt)'\ejm)*M.muj(n,jj,Gtt)'*dtdDphinjm;
        %                     end
        %                     
        %                     % debug
        %                     dtdmunj = dy(mujoffset+dim*(jj-1)+(1:dim),n);
        %                     dtdDphin = reshape(dy(Dphioffset+(1:dim*dim),n),dim,dim);
        %                     dtdmunjCheck = inv(M.Dphi(n,Gtt))'*dDphin'*inv(M.Dphi(n,Gtt))'*M.dtDphi(n,Gtt)'*M.muj(n,jj,Gtt)...
        %                         - inv(M.Dphi(n,Gtt))'*dtdDphin'*M.muj(n,jj,Gtt)...
        %                         - inv(M.Dphi(n,Gtt))'*M.dtDphi(n,Gtt)'*dmunj;
        %                     e = norm(dtdmunj-dtdmunjCheck);
        %                     if e > 1e-6
        %                         fprintf('dtdmunj error (t=%f): %f\n',t,e);
        %                     end                    
        %                 end            
        %             end
        %         end
        %     end
        % 
        %     dy = reshape(dy,gradCSP*L,1);
        % end
        % 
        % % debug (forward systems)
        % % muj
        % 
        % function dy = Gmuj(t,ytt)
        %     ytt = reshape(ytt,dim*dim,L);
        %     
        %     dy = zeros(size(ytt));
        %     
        %     Gtt = reshape(deval(Gt,t),dim*(1+dim+1),L);
        %     
        %     for n = 1:L % update all particles        
        %         for jj = 1:dim
        %             munj = ytt(dim*(jj-1)+(1:dim),n);
        %             a0nj = rhoj0(dim*(jj-1)+(1:dim),n);
        %             munjCheck = M.Dphi(n,Gtt)'\a0nj;
        %             e = norm(munj-munjCheck);
        %             if e > 1e-5
        %                 fprintf('munj error (t=%f): %f\n',t,e);
        %             end   
        %             assert(e < 1e-5);
        % 
        % %             dy(dim*(jj-1)+(1:dim),n) = dy(dim*(jj-1)+(1:dim),n) ...
        % %                 + M.Dphi(n,Gtt)'\(M.dtDphi(n,Gtt)'*munj);
        %             dy(dim*(jj-1)+(1:dim),n) = dy(dim*(jj-1)+(1:dim),n) ...
        %                 - (inv(M.Dphi(n,Gtt))*M.dtDphi(n,Gtt)*inv(M.Dphi(n,Gtt)))'*a0nj;
        %         end
        %     end
        % 
        %     dy = reshape(dy,dim*dim*L,1);
        % end
        % options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        % wtmuj = ode45(@Gmuj,[0 1],reshape(rhoj0,dim*dim*L,1),options);
        % % check full forward system
        % function doForward(da0,da0j)
        %     initial = reshape([zeros(dim*(1+dim),L); da0; da0j;],gradCSP*L,1);
        %     
        % %     xx = load('debug.mat');
        % %     q1 = Gforward1(xx.t,xx.ytt);
        % %     q2 = Gforward2(xx.t,xx.ytt);
        % %     reshape(q1-q2,[],L)    
        % 
        %     options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        %     wt12 = ode45(@Gforward,[0 1],initial,options);
        %     wt2 = ode45(@Gforward2,[0 1],initial,options);
        %     wt1 = ode45(@Gforward1,[0 1],initial,options);            
        % 
        %     norm(deval(wt1,1)-deval(wt2,1))
        %     %assert(norm(deval(wt1,1)-deval(wt2,1)) < 1e-5);    
        % end
        % for dn = 1:L
        %     for dj = 1:dim
        %         da0 = zeros(dim,L);
        %         da0(dj,dn) = 1;
        %         
        %         % da0
        %         da0j = zeros(dim*dim,L);
        %         
        %         doForward(da0,da0j);
        % 
        %         % da0j
        %         da0 = zeros(dim,L);       
        %         for djm = 1:dim
        %             da0j = zeros(dim*dim,L);
        %             da0j(dim*(dj-1)+djm,dn) = 1;
        %             
        %             doForward(da0,da0j);
        %         end
        %     end
        % end

        % integrate
        options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        if cdim == dim
            w1 = [reshape(v1,CSP,L); zeros(dim^2,L)];
        else
            assert(dim == 2 && cdim == 3); % shift from 2d to 3d
            v1 = reshape(v1,CSP,L);
            w1 = zeros(cgradCSP,L);
            w1(1:dim,:) = v1(1:dim,:);
            w1(cdim+(1:cdim:cdim*dim),:) = v1(dim+(1:dim:dim^2),:);
            w1(cdim+(2:cdim:cdim*dim),:) = v1(dim+(2:dim:dim^2),:);
            w1(cdim+cdim^2+(1:dim),:) = v1(dim+dim^2+(1:dim),:);
            v1 = reshape(v1,CSP*L,1);
            % rest is zero
        end
        wt = ode45(@Gc,[0 1],w1,options); % solve backwards, fast native
%         wt = ode45(@G,[0 1],w1,options); % solve backwards, slooow matlab
        assert(wt.x(end) == 1);
        w0 = reshape(deval(wt,1),cgradCSP,L);
        if dim == cdim
            v0 = [w0(1:dim,:); w0(muoffset+1:end,:)]; % phi, mu, and muj parts
        else
            assert(dim == 2 && cdim == 3); % shift from 3d to 2d
            v0 = zeros(CSP,L);
            v0(1:dim,:) = w0(1:dim,:);
            v0(dim+(1:dim),:) = w0(muoffset+(1:dim),:);
            v0(2*dim+(1:dim:dim*dim),:) = w0(mujoffset+(1:cdim:cdim*dim),:);
            v0(2*dim+(2:dim:dim*dim),:) = w0(mujoffset+(2:cdim:cdim*dim),:);
        end                 
        v0 = reshape(v0,CSP*L,1);
    end

gradTransport = @lgradTransport;

end