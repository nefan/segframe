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

function pathEnergy = getPointPathEnergy(lddmmoptions)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim;
L = lddmmoptions.L;
R = lddmmoptions.R;
CSP = lddmmoptions.CSP;
cCSP = lddmmoptions.cCSP;
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;
order = lddmmoptions.order;

ks = dkernelsGaussian(cdim);

    function [E,v0] = lpathEnergy(x,Gtt,varargin)
        tend = 1;
        if size(varargin,2) > 0
            tend = varargin{1};
        end
        
        function Et = G(tt)
            t = intTime(tt,false,lddmmoptions);
            Gt = reshape(deval(Gtt,t),cCSP,L);

            Et = 0;

            q0 = tensor(Gt(1:cdim,:),[cdim L],'ai');
            switch order
                case 0
                    mu0 = tensor(Gt((cdim+1):(cdim+cdim),:),[cdim L],'ai');            
                case 1
                    q1 = tensor(Gt((cdim+1):(cdim+cdim^2),:),[cdim cdim L],'abi');
                    mu0 = tensor(Gt((cdim+cdim^2+1):(2*cdim+cdim^2),:),[cdim L],'ai');
                    mu1 = tensor(Gt((2*cdim+cdim^2+1):(2*cdim+2*cdim^2),:),[cdim cdim L],'abi');
                case 2
                    q1 = tensor(Gt((cdim+1):(cdim+cdim^2),:),[cdim cdim L],'abi');
                    q2 = tensor(Gt((cdim+cdim^2+1):(cdim+cdim^2+cdim^3),:),[cdim cdim cdim L],'abgi');
                    mu0 = tensor(Gt((cdim+cdim^2+cdim^3+1):(2*cdim+cdim^2+cdim^3),:),[cdim L],'ai');
                    mu1 = tensor(Gt((2*cdim+cdim^2+cdim^3+1):(2*cdim+2*cdim^2+cdim^3),:),[cdim cdim L],'abi');
                    mu2 = tensor(Gt((2*cdim+2*cdim^2+cdim^3+1):(2*cdim+2*cdim^2+2*cdim^3),:),[cdim cdim cdim L],'abgi');
            end

            % compute kernel and derivatives
            switch order
                case 0
                    [Ks] = ks.TKs(q0,q0,scales,scaleweight);
                case 1
                    [Ks,D1Ks,D2Ks] = ks.TKs(q0,q0,scales,scaleweight);
                case 2
                    [Ks,D1Ks,D2Ks] = ks.TKs(q0,q0,scales,scaleweight);
            end        

            % q
            Et = Et + sum(ttov(tcntr(tproddiag(tproddiag(tind(mu0,'aj'),Ks,'j'),mu0,'i'),'a')));
            if order >= 1
                Et = Et - sum(ttov(tcntr(tproddiag(tcntr(tproddiag(tind(mu1,'abj'),D2Ks,'j'),'b'),mu1,'i'),'a')));
                Et = Et + 2*sum(ttov(tcntr(tproddiag(tcntr(tproddiag(tind(mu1,'abj'),D1Ks,'j'),'b'),mu0,'i'),'a')));
            end
            if order >= 2
                assert(false);
            end
            Et

        end

        E = integrate(@G,0,tend);
        assert(E >= 0);
        
        v0 = [];
    end

switch order
    case 0
        pathEnergy = getPointPathEnergyOrder0(lddmmoptions);
    case 1
        pathEnergy = getPointPathEnergyOrder1(lddmmoptions);
%         pathEnergy = @lpathEnergy;
    case 2
        pathEnergy = @lpathEnergy;
end


end