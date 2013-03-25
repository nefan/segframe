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

function pathEnergy = getPointPathEnergyOrder1(lddmmoptions)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim;
L = lddmmoptions.L;
R = lddmmoptions.R;
CSP = lddmmoptions.CSP;
cCSP = lddmmoptions.cCSP;
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;

ks = dkernelsGaussian(cdim);

    function [E v0] = lpathEnergy(x,Gtt)
        
        x = reshape(x,(1+dim)*dim,L);
        rhoj0 = x(dim+(1:dim^2),:);
        if dim ~= cdim
            rhoj0 = reshape(rhoj2dTo3dOrder1(rhoj0,lddmmoptions),cdim^2,L);
        end      
        
        
        function v = gradE()    
            assert(R == 1);
            s = 1;
            v = zeros(cdim*(1+cdim),L);
            G0 = reshape(deval(Gtt,0),cdim*(1+cdim+1),L);
            for n = 1:L
                xn = G0(1:cdim,n);
                an = G0(cdim*(1+cdim)+(1:cdim),n); 

                for ll = 1:L % particle
                    xl = G0(1:cdim,ll);
                    al = G0(cdim*(1+cdim)+(1:cdim),ll);  
                    Ks = ks.Ks(xl,xn,scales(s),scaleweight(s));
                    N1 = ks.N1Ks(xl,xn,scales(s),scaleweight(s));
                    D2N1 = ks.D2N1Ks(xl,xn,scales(s),scaleweight(s));

                    % an
                    v(1:cdim,n) = v(1:cdim,n) + 2*Ks*al;

                    for j = 1:cdim
                        alj = rhoj0(cdim*(j-1)+(1:cdim),ll);             
                        % an
                        v(1:cdim,n) = v(1:cdim,n) + 2*N1(j)*alj;

                        % anj
                        v(cdim+cdim*(j-1)+(1:cdim),n) = v(cdim+cdim*(j-1)+(1:cdim),n) - 2*N1(j)*an;

                        for jm = 1:cdim
                            aljm = rhoj0(cdim*(j-1)+(1:cdim),ll);             

                            % anj
                            v(cdim+cdim*(j-1)+(1:cdim),n) = v(cdim+cdim*(j-1)+(1:cdim),n) + 2*D2N1(jm,j)*aljm;
                        end
                    end            
                end
            end

            v = reshape(v,cdim*(1+cdim)*L,1);
        end
        
        function Et = G(tt)
            t = intTime(tt,false,lddmmoptions);
            Gt = reshape(deval(Gtt,t),cCSP,L);

            Et = 0;

            for s = 1:R
                for i = 1:L
                    xi = Gt(1:cdim,i);
                    Dphii = reshape(Gt(cdim+(1:cdim*cdim),i),cdim,cdim);
                    mui = Gt(cdim+cdim^2+(1:cdim),i);
                    ai = mui;

                    for l = 1:L
                        xl = Gt(1:cdim,l);
                        Dphil = reshape(Gt(cdim+(1:cdim^2),l),cdim,cdim);
                        mul = Gt(cdim+cdim^2+(1:cdim),l);
                        al = mul;

                        Et = Et + ai'*ks.Ks(xl,xi,scales(s),scaleweight(s))*al;

                        N1Kli = ks.N1Ks(xl,xi,scales(s),scaleweight(s));                
                        N1Kil = ks.N1Ks(xi,xl,scales(s),scaleweight(s));                
                        D2N1 = ks.D2N1Ks(xl,xi,scales(s),scaleweight(s));
                        assert(norm(D2N1-D2N1') < 1e-10);
                        assert(norm(D2N1-ks.D1N2Ks(xl,xi,scales(s),scaleweight(s))) < 1e-10);
                        assert(norm(D2N1-ks.D2N1Ks(xi,xl,scales(s),scaleweight(s))) < 1e-10);

                        for j = 1:cdim
                            alj = zeros(cdim,1);
                            aij = zeros(cdim,1);
                            for k=1:cdim
                                mulk = Dphil'\rhoj0(cdim*(k-1)+(1:cdim),l);
                                alj = alj+Dphil(j,k)*mulk;

                                muik = Dphii'\rhoj0(cdim*(k-1)+(1:cdim),i);
                                aij = aij+Dphii(j,k)*muik;                        
                            end

                            Et = Et + N1Kli(j)*ai'*alj+N1Kil(j)*al'*aij;
                            %Et = Et + 2*N1Kli(j)*ai'*alj;
                            assert(N1Kli(j) == -N1Kil(j));

                            for jm = 1:cdim
                                aijm = zeros(cdim,1);
                                for k=1:cdim
                                    muik = Dphii'\rhoj0(cdim*(k-1)+(1:cdim),i);
                                    aijm = aijm+Dphii(jm,k)*muik;
                                end

                                Et = Et + D2N1(j,jm)*aijm'*alj;                        
                            end
                        end
                    end
                end
            end

        end

        % Epath = integrate(@Gc,0,1); % fast C version
        E = integrate(@G,0,1); % sloow matlab version
        assert(E >= 0);
        
        % warning: x-parts not included here!
        w0 = reshape(gradE(),cdim^2+cdim,L);
        if dim == cdim
            v0 = [zeros(dim,L); w0];
        else
            assert(dim == 2 && cdim == 3); % shift from 3d to 2d
            v0 = zeros(CSP,L);
            v0(dim+(1:dim),:) = w0(1:dim,:);
            v0(2*dim+(1:dim:dim*dim),:) = w0(cdim+(1:cdim:cdim*dim),:);
            v0(2*dim+(2:dim:dim*dim),:) = w0(cdim+(2:cdim:cdim*dim),:);            
        end
        v0 = reshape(v0,CSP*L,1);
    end

pathEnergy = @lpathEnergy;

end