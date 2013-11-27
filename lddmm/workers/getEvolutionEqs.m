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

function evoG = getEvolutionEqs(lddmmoptions,varargin)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim; % computations performed in cdim
L = lddmmoptions.L;
R = lddmmoptions.R;
order = lddmmoptions.order;
cCSP = lddmmoptions.cCSP;
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;
epsilon = lddmmoptions.epsilon;

[Ks D1Ks D2Ks] = gaussianKernels(); % depreceated
ks = dkernelsGaussian(cdim);

assert(order >= 0 && order <= 2); % supported orders

rhoj0 = [];
if size(varargin,2) > 0
    rhoj0 = varargin{1};
end

    function drho = G0M(tt,rhot) % Matrix version
        t = intTime(tt,false,lddmmoptions);
        
        rhots = reshape(rhot,cCSP,L);
        xt = rhots(1:cdim,:);
        at = rhots((cdim+1):end,:);
        
        function v = Mf(pmi,pml,i,l)
            v = zeros(cdim);
            sl = pml-1;
            
            if pmi == 1 % position, row
               if pml == 1 % position, col
               else % momentum, col
                   sl = pml-1;
                   v = ks.Ks(xt(:,i),xt(:,l),scales(sl),scaleweight(sl))*eye(cdim);
               end
            else % momentum, row
               if pml == 1 % position, col
               else % momentum, col                   
                   v = -ks.N1Ks(xt(:,i),xt(:,l),scales(sl),scaleweight(sl))*at(:,i)';
               end                
            end
            
            v = {mtov(v)};   
        end
        
        M = reshape(cell2mat(mmap(@Mf,1+R,1+R,L,L)),cdim,cdim,1+R,1+R,L,L);
        M = reshape(permute(M,[1 3 5 2 4 6]),cdim*(1+R)*L,cdim*(1+R)*L);
        
        drho = M*rhot;
        
        drho = intResult(drho,false,lddmmoptions);

        % debug
        if getOption(lddmmoptions,'testC')
            drho2 = G0(tt,rhot);
            if norm(drho-drho2) > 10e-12
                1;
            end
            assert(norm(drho-drho2) < 10e-12);
        end
    end

    function xt = Gindex(tt,xx) % tensor version
        t = intTime(tt,false,lddmmoptions);
        
        x = reshape(xx,cCSP,L);
        q0 = tensor(x(1:cdim,:),[cdim L],'ai');
        switch order
            case 0
                mu0 = tensor(x((cdim+1):(cdim+cdim),:),[cdim L],'ai');            
            case 1
                q1 = tensor(x((cdim+1):(cdim+cdim^2),:),[cdim cdim L],'abi');
                mu0 = tensor(x((cdim+cdim^2+1):(2*cdim+cdim^2),:),[cdim L],'ai');
                mu1 = tensor(x((2*cdim+cdim^2+1):(2*cdim+2*cdim^2),:),[cdim cdim L],'abi');
            case 2
                assert(false)
        end
        
        % compute kernel and derivatives
        switch lddmmoptions.order
            case 0
                [Ks,D1Ks] = ks.TKs(q0,scales,scaleweight);
            case 1
                [Ks,D1Ks,D2Ks,D3Ks] = ks.TKs(q0,scales,scaleweight);
            case 2
                assert(false);
        end        

        % chi 
        e1 = tprodcntr(mu0,2,D1Ks,2);
        if lddmmoptions.order == 1 
            e1 = tsum(e1,tcntr(tprodcntr(mu1,2,D2Ks,4),2,4));
            e2 = tsum(tprodcntr(mu0,2,D2Ks,2),tcntr(tprodcntr(mu1,2,D3Ks,5),2,4));
        end        
        if lddmmoptions.order == 2             
        end                
        
        % mu
        mu0t = tscalar(-1,tshift(tcntr(tproddiag(mu0,2,e1,2),1,3),[2 1]));
        if lddmmoptions.order == 1 
            mu0t = tsum(mu0t,tshift(tcntr(tcntr(tproddiag(mu1,3,e2,2),2,6),1,3),[2 1]));
            mu0t.indices = 'ai';

            T = tproddiag(mu1,3,e1,2);
            mu1t = tsub(tshift(tcntr(T,2,5),[1 3 2]),tshift(tcntr(T,1,4),[3 1 2]));
            mu1t.indices = 'abi';
        end
        if lddmmoptions.order == 2             
        end                

        % q
        q0t = tprodcntr(mu0,2,Ks,2);
        if lddmmoptions.order == 1 
            q0t = tsum(q0t,tcntr(tprodcntr(mu1,2,D1Ks,3),2,4));

            q1t = tshift(tdiag(tprodcntr(e1,3,q1,1),2,4),[1 3 2]);
        end                
        if lddmmoptions.order == 2             
        end                
        
        switch order
            case 0
                xt = [q0t.T; mu0t.T];
            case 1
                xt = [q0t.T; reshape(q1t.T,cdim^2,L); mu0t.T; reshape(mu1t.T,cdim^2,L)];
            case 2
                assert(false);
        end
        
        xt = reshape(xt,[],1);
        
        xt = intResult(xt,false,lddmmoptions);


%         % debug
%         if getOption(lddmmoptions,'testC')
%             if lddmmoptions.order == 0
%                 drho2 = G0(tt,xx);
%             else if lddmmoptions.order == 1
%                 drho2 = G1(tt,xx);
%             else if lddmmoptions.order == 2
%                     assert(false);
%             end
%             end
%             end
% 
%             xtdebug = reshape(xt,[],L);
%             xtdebug = reshape(xtdebug(1:(cdim+cdim^2+cdim),:),[],1);
%             drho2 = reshape(drho2,[],L);
%             drho2 = reshape(drho2(1:(cdim+cdim^2+cdim),:),[],1);
%             
%             if norm(xtdebug-drho2) > 1e-7
%                 norm(xtdebug-drho2)
%                 tt
%                 1;                
%             end
%             assert(norm(xtdebug-drho2) < 1e-12);
%         end
    end

    function drho = G0c(tt,rhot) % wrapper for C version
        t = intTime(tt,false,lddmmoptions);
        
        drho = fastPointPathOrder0(t,rhot,L,R,cdim,scales.^2,scaleweight.^2);
        
        drho = intResult(drho,false,lddmmoptions);

        % debug
        if getOption(lddmmoptions,'testC')
            drho2 = G0(tt,rhot);
            if norm(drho-drho2) > 10e-12
                1;
            end
            assert(norm(drho-drho2) < 10e-12);
        end
    end

    function drho = G0(tt,rhot)  % slooow version
        rhot = reshape(rhot,cCSP,L);
        t = intTime(tt,false,lddmmoptions);

        drho = zeros(size(rhot));

        for i = 1:L % particle
            xi = rhot(1:cdim,i);

            for l = 1:L % particle
                xl = rhot(1:cdim,l);

                ximxl = xi-xl;

                for sl = 1:R % scale
                    d1ksl = D1Ks(xi,xl,scales(sl),scaleweight(sl));
                    TWOd1kslximxl = 2*d1ksl*ximxl;

                    % position
                    drho(1:cdim,i) = drho(1:cdim,i)+Ks(xi,xl,scales(sl),scaleweight(sl))*rhot(cdim*(sl-1)+(1+cdim:2*cdim),l);

                    for si = 1:R % scale
                        % momentum
                        drho(cdim*(si-1)+(1+cdim:2*cdim),i) = drho(cdim*(si-1)+(1+cdim:2*cdim),i) ...
                            - TWOd1kslximxl*rhot(cdim*(sl-1)+(1+cdim:2*cdim),l)'*rhot(cdim*(si-1)+(1+cdim:2*cdim),i);
                    end                
                end
            end
        end

        drho = intResult(drho,false,lddmmoptions);
        drho = reshape(drho,L*cCSP,1);
    end

        function dGt = G1c(tt,Gt) % wrapper for cpu version of G
            t = intTime(tt,false,lddmmoptions);

            dGt = fastPointPathOrder1(t,Gt,reshape(rhoj0,[],1),L,R,cdim,scales.^2,scaleweight.^2);
            
            dGt = intResult(dGt,false,lddmmoptions);

            % debug
            if getOption(lddmmoptions,'testC')
                dGt2 = G1(tt,Gt);  
                assert(norm(dGt-dGt2) < epsilon);
            end
        end           
        
%         function drho = G1M(tt,Gt) % Matrix version
%             t = intTime(tt,false,lddmmoptions);
% 
%             rhots = reshape(rhot,cCSP,L);
%             xt = rhots(1:cdim,:);
%             at = rhots((cdim+1):end,:);
% 
%             function v = Mf(pmi,pml,i,l)
%                 v = zeros(cdim);
% 
%                 if pmi == 1 % position, row
%                    if pml == 1 % position, col
%                    else % momentum, col
%                        sl = pml-1;
%                        v = Ks(xt(:,i),xt(:,l),scales(sl),scaleweight(sl))*eye(cdim);
%                    end
%                 else % momentum, row
%                    if pml == 1 % position, col
%                    else % momentum, col
%                        sl = pml-1;
%                        v = -2*D1Ks(xt(:,i),xt(:,l),scales(sl),scaleweight(sl))*(xt(:,i)-xt(:,l)) ...
%                              *at(:,i)';
%                    end                
%                 end
% 
%                 v = {v};   
%             end
% 
%             M = reshape(cell2mat(mmap(@Mf,1+R,L,1+R,L)),cdim*(1+R),cdim*(1+R),L,L);
%             M = reshape(permute(M,[1 3 2 4]),cdim*(1+R)*L,cdim*(1+R)*L);
% 
%             drho = M*rhot;
% 
%             drho = intResult(drho,false,lddmmoptions);
% 
%             % debug
%             if getOption(lddmmoptions,'testC')
%                 drho2 = G1(tt,rhot);
%                 if norm(drho-drho2) > 10e-12
%                     1;
%                 end
%                 assert(norm(drho-drho2) < 10e-12);
%             end
%         end        
    
        function dGt = G1(tt,Gt)  % slooow version
            Gt = reshape(Gt,cCSP,L);
            t = intTime(tt,false,lddmmoptions);

            dGt = zeros(size(Gt));

            for n = 1:L % particle
                xn = Gt(1:cdim,n);
                Dphin = reshape(Gt(cdim+(1:cdim^2),n),cdim,cdim);

                for i = 1:L % particle
                    xi = Gt(1:cdim,i);
                    Dphii = reshape(Gt(cdim+(1:cdim^2),i),cdim,cdim);            

                    for si = 1:R % scale
                        mui = Gt(cdim+cdim^2+(1:cdim),i);

                        % position
                        dGt(1:cdim,n) = dGt(1:cdim,n)+...
                            ks.Ks(xi,xn,scales(si),scaleweight(si))*mui;
                        for j = 1:cdim
                            muij = Dphii'\rhoj0(cdim*(j-1)+(1:cdim),i);
                            dGt(1:cdim,n) = dGt(1:cdim,n)...
                                +ks.N1Ks(xi,xn,scales(si),scaleweight(si))'*Dphii(:,j)*muij;
                        end
                        % dPhi
                        for c = 1:cdim
                            dGt(cdim+cdim*(c-1)+(1:cdim),n) = dGt(cdim+cdim*(c-1)+(1:cdim),n)...
                                +ks.N2Ks(xi,xn,scales(si),scaleweight(si))'*Dphin(:,c)*mui;
                            for j = 1:cdim
                                muij = Dphii'\rhoj0(cdim*(j-1)+(1:cdim),i);
                                dGt(cdim+cdim*(c-1)+(1:cdim),n) = dGt(cdim+cdim*(c-1)+(1:cdim),n)...
                                    +(ks.D1N2Ks(xi,xn,scales(si),scaleweight(si))*Dphii(:,j))'*Dphin(:,c)*muij;
                            end                
                        end

                        for sn = 1:R % scale
                            mun = Gt(cdim+cdim^2+(1:cdim),n);

                            % mu
                            dGt(cdim+cdim^2+(1:cdim),n) = dGt(cdim+cdim^2+(1:cdim),n)...
                                -mun'*mui*ks.N2Ks(xi,xn,scales(si),scaleweight(si));
                            for j = 1:cdim
                                muij = Dphii'\rhoj0(cdim*(j-1)+(1:cdim),i);
                                munj = Dphin'\rhoj0(cdim*(j-1)+(1:cdim),n);
                                dGt(cdim+cdim^2+(1:cdim),n) = dGt(cdim+cdim^2+(1:cdim),n)...
                                    -munj'*mui*ks.D2N2Ks(xi,xn,scales(si),scaleweight(si))*Dphin(:,j) ...
                                    -mun'*muij*ks.D1N2Ks(xi,xn,scales(si),scaleweight(si))*Dphii(:,j);
                                    %-(munj'*mui-mun'*muij)*D2N2Ks(xi,xn,scales(si),scaleweight(si))*Dphin(:,j);                       
                                for jm = 1:cdim
                                    munjm = Dphin'\rhoj0(cdim*(jm-1)+(1:cdim),n);
                                    dGt(cdim+cdim^2+(1:cdim),n) = dGt(cdim+cdim^2+(1:cdim),n)...
                                        -munjm'*muij*ks.D2D1N2Ksa(Dphii(:,j),xi,xn,scales(si),scaleweight(si))*Dphin(:,jm);
                                end                                                            
                            end                                    
                        end                
                    end
                end

            end

            dGt = intResult(dGt,false,lddmmoptions);            
            dGt = reshape(dGt,L*cCSP,1);
        end

switch order
    case 0
        evoG = @Gindex;
    case 1
        evoG = @Gindex;        
    case 2
        assert(false);
end

end