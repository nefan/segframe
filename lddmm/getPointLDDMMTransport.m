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

function transport = getPointLDDMMTransport(pointPath,moving,lddmmoptions)

dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim; % computations performed in cdim
L = lddmmoptions.L;
R = lddmmoptions.R;
order = lddmmoptions.order;
CSP = lddmmoptions.CSP;
cCSP = lddmmoptions.cCSP;
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;

ks = dkernelsGaussian(cdim);

    function [grid1,gridt] = lpathtransport(grid0,Gtt,backwards,tend)
        Ngrid = numel(grid0{1});
        
        function gridt = Gindex(t,grid)
            t = intTime(t,backwards,lddmmoptions);    

            grid = tensor(grid,[cdim Ngrid],'ai');

            Gt = reshape(deval(Gtt,t),cCSP,L);

            q0 = tensor(Gt(1:cdim,:),[cdim L],'ai');
            switch order
                case 0
                    mu0 = tensor(Gt((cdim+1):(cdim+cdim),:),[cdim L],'ai');            
                case 1
                    mu0 = tensor(Gt((cdim+cdim^2+1):(2*cdim+cdim^2),:),[cdim L],'ai');
                    mu1 = tensor(Gt((2*cdim+cdim^2+1):(2*cdim+2*cdim^2),:),[cdim cdim L],'abi');
                case 2
                    mu0 = tensor(Gt((cdim+cdim^2+cdim^3+1):(2*cdim+cdim^2+cdim^3),:),[cdim L],'ai');
                    mu1 = tensor(Gt((2*cdim+cdim^2+cdim^3+1):(2*cdim+2*cdim^2+cdim^3),:),[cdim cdim L],'abi');
                    mu2 = tensor(Gt((2*cdim+2*cdim^2+cdim^3+1):(2*cdim+2*cdim^2+2*cdim^3),:),[cdim cdim cdim L],'abgi');
            end

            % compute kernel and derivatives
            switch order
                case 0
                    [Ks] = ks.TKs(grid,q0,scales,scaleweight);
                case 1
                    [Ks,D1Ks] = ks.TKs(grid,q0,scales,scaleweight);
                case 2
                    [Ks,D1Ks,D2Ks,D3Ks] = ks.TKs(grid,q0,scales,scaleweight);
            end                     

            % gridt (t derivative)
            gridt = tprodcntr(tind(mu0,'aj'),Ks,'j');
            if order >= 1 
                gridt = tsub(gridt,tcntr(tprodcntr(tind(mu1,'abj'),D1Ks,'j'),'b'));
            end                
            if order >= 2
                gridt = tsum(gridt,tcntr(tcntr(tprodcntr(tind(mu2,'abgj'),D2Ks,'j'),'g'),'b'));
            end
            
            gridt = ttov(gridt);                             
            gridt = intResult(gridt,backwards,lddmmoptions);

        end
        
        siz = size(grid0{1});
        g0 = reshape([reshape(grid0{1},1,Ngrid); reshape(grid0{2},1,Ngrid); reshape(grid0{3},1,Ngrid)],3*Ngrid,1);
        gridt = ode45(@Gindex,[0 tend],g0);
        g1 = reshape(deval(gridt,tend),3,Ngrid);
        grid1 = cell(1);
        grid1{1} = reshape(g1(1,:),siz);
        grid1{2} = reshape(g1(2,:),siz);
        grid1{3} = reshape(g1(3,:),siz);
    end

    function [res,Gt] = ltransport(x, points, varargin)
        backwards = false;
        if size(varargin,2) > 0
            backwards = varargin{1};
        end
        tend = 1;
        if size(varargin,2) > 1
            tend = varargin{2};
        end
        
        ic = pathIC(x,moving,lddmmoptions);        
        Gt = pointPath(ic,tend); % rhot is a solver structure
        
        if iscell(points)
            g0 = points;
        else
            g0{1} = points(1,:);
            g0{2} = points(2,:);
            if dim == 2
                g0{3} = zeros(size(g0{1}));
            else
                g0{3} = points(3,:);
            end
        end
        
        switch order
            case 0
%                 g1 = shootgrid(g0,Gt,lddmmoptions,backwards,tend);
                g1 = lpathtransport(g0,Gt,backwards,tend);
            case 1
%                 g1 = shootgridDKernels(x,g0,Gt,lddmmoptions,backwards);
                g1 = lpathtransport(g0,Gt,backwards,tend);
            case 2
                g1 = lpathtransport(g0,Gt,backwards,tend);
        end
        if iscell(points)
            res = g1;
        else
            switch dim
                case 2
                    res = [g1{1}; g1{2}];
                case 3
                    res = [g1{1}; g1{2}; g1{3}];
            end
        end
    end

transport = @ltransport;

end
