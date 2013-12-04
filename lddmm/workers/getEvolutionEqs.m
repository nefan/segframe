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

ks = dkernelsGaussian(cdim);

assert(order >= 0 && order <= 2); % supported orders

rhoj0 = []; % depreceated
if size(varargin,2) > 0
    rhoj0 = varargin{1};
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
                q1 = tensor(x((cdim+1):(cdim+cdim^2),:),[cdim cdim L],'abi');
                q2 = tensor(x((cdim+cdim^2+1):(cdim+cdim^2+cdim^3),:),[cdim cdim cdim L],'abgi');
                mu0 = tensor(x((cdim+cdim^2+cdim^3+1):(2*cdim+cdim^2+cdim^3),:),[cdim L],'ai');
                mu1 = tensor(x((2*cdim+cdim^2+cdim^3+1):(2*cdim+2*cdim^2+cdim^3),:),[cdim cdim L],'abi');
                mu2 = tensor(x((2*cdim+2*cdim^2+cdim^3+1):(2*cdim+2*cdim^2+2*cdim^3),:),[cdim cdim cdim L],'abgi');
        end
        
        % compute kernel and derivatives
        switch order
            case 0
                [Ks,D1Ks] = ks.TKs(q0,q0,scales,scaleweight);
            case 1
                [Ks,D1Ks,D2Ks,D3Ks] = ks.TKs(q0,q0,scales,scaleweight);
            case 2
                [Ks,D1Ks,D2Ks,D3Ks,D4Ks,D5Ks] = ks.TKs(q0,q0,scales,scaleweight);
        end        

        % chi 
        e1 = tprodcntr(tind(mu0,'aj'),D1Ks,'j');
        if order >= 1 
            e1 = tsum(e1,tcntr(tprodcntr(tind(mu1,'agj'),D2Ks,'g'),'j'));
            e2 = tsum(tprodcntr(tind(mu0,'aj'),D2Ks,'j'),tcntr(tprodcntr(tind(mu1,'adj'),D3Ks,'d'),'j'));
        end        
        if order >= 2
            e1 = tsum(e1,tcntr(tcntr(tprodcntr(tind(mu2,'agdj'),D3Ks,'d'),'g'),'j'));
            e2 = tsum(e2,tcntr(tcntr(tprodcntr(tind(mu2,'adej'),D4Ks,'e'),'d'),'j'));            
            e3 = tsum(tsum(tprodcntr(tind(mu0,'aj'),D3Ks,'j'),tcntr(tprodcntr(tind(mu1,'aej'),D4Ks,'e'),'j')),...
                      tcntr(tcntr(tprodcntr(tind(mu2,'aepj'),D5Ks,'p'),'e'),'j'));
                  e2.T = 0*e2.T;
%                   e3.T = 0*e3.T;
        end                
        
        % mu
        mu0t = tscalar(-1,tshift(tcntr(tproddiag(tind(mu0,'bi'),tind(e1,'bia'),'i'),'b'),[2 1]));
        if order >= 1 
            mu0t = tsum(mu0t,tshift(tcntr(tcntr(tproddiag(tind(mu1,'bgi'),tind(e2,'biga'),'i'),'g'),'b'),[2 1]));

            T = tproddiag(mu1,e1,'i');
            mu1t = tsub(tshift(tcntr(tind(T,'agibg'),'g'),[1 3 2]),tshift(tcntr(tind(T,'gbiga'),'g'),[3 1 2]));
        end
        if order >= 2             
            mu0t = tsub(mu0t,tshift(tcntr(tcntr(tcntr(tproddiag(tind(mu2,'bgdi'),tind(e3,'bigda'),'i'),'d'),'g'),'b'),[2 1]));
            
            T = tproddiag(mu2,e2,'i');
            mu1t = tsub(tsub(tsum(mu1t,tshift(tcntr(tcntr(tind(T,'adgibdg'),'g'),'d'),[1 3 2])),...
                             tshift(tcntr(tcntr(tind(T,'dbgidag'),'g'),'d'),[3 1 2])),...
                        tshift(tcntr(tcntr(tind(T,'dgbidga'),'g'),'d'),[3 1 2]));
            
            T = tproddiag(mu2,e1,'i');
            mu2t = tshift(tsub(tsum(tcntr(tind(T,'adgibd'),'d'),...
                                    tcntr(tind(T,'agdibd'),'d')),...
                               tcntr(tind(T,'dbgida'),'d')),[4 1 2 3]);
        end                

        % q
        q0t = tprodcntr(tind(mu0,'aj'),Ks,'j');
        if order >= 1 
            q0t = tsum(q0t,tcntr(tprodcntr(tind(mu1,'abj'),D1Ks,'j'),'b'));
            q1t = tshift(tcntr(tproddiag(tind(e1,'aig'),tind(q1,'gbi'),'i'),'g'),[1 3 2]);
        end                
        if order >= 2
            q0t = tsum(q0t,tcntr(tcntr(tprodcntr(tind(mu2,'abgj'),D2Ks,'j'),'g'),'b'));
            q2t = tshift(tsum(tcntr(tproddiag(tcntr(tproddiag(tind(e2,'aide'),tind(q1,'dbi'),'i'),'d'),tind(q1,'egi'),'i'),'e'),...
                              tcntr(tproddiag(tind(e1,'aid'),tind(q2,'dbgi'),'i'),'d')),[1 3 4 2]);
        end                
        
        switch order
            case 0
                xt = [q0t.T; mu0t.T];
            case 1
                xt = [q0t.T; reshape(q1t.T,cdim^2,L); mu0t.T; reshape(mu1t.T,cdim^2,L)];
            case 2
                xt = [q0t.T; reshape(q1t.T,cdim^2,L); reshape(q2t.T,cdim^3,L); mu0t.T; reshape(mu1t.T,cdim^2,L); reshape(mu2t.T,cdim^3,L)];
        end
        
        xt = reshape(xt,[],1);
        
        xt = intResult(xt,false,lddmmoptions);


%         % debug
%         if getOption(lddmmoptions,'testC')
%             if order == 0
%                 drho2 = G0(tt,xx);
%             else if order == 1
%                 drho2 = G1(tt,xx);
%             else if order == 2
%                     assert(false);
%             end
%             end
%             end
% 
%             xtdebug = xt;
% %             xtdebug = reshape(xt,[],L);
% %             xtdebug = reshape(xtdebug(1:(cdim+cdim^2+cdim),:),[],1);
% %             drho2 = reshape(drho2,[],L);
% %             drho2 = reshape(drho2(1:(cdim+cdim^2+cdim),:),[],1);
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

        [Ks D1Ks D2Ks] = gaussianKernels();

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
        evoG = @Gindex;        
end

end