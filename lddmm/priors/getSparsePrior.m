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

function prior = getSparsePrior(lddmmoptions)

dim = lddmmoptions.dim;
L = lddmmoptions.L;
R = lddmmoptions.R;
CSP = lddmmoptions.CSP;
order = lddmmoptions.order;
scales = lddmmoptions.scales;

alpha = lddmmoptions.sparseoptions.alpha;
c = lddmmoptions.sparseoptions.c;
lambdasp = lddmmoptions.sparseoptions.lambdasp;

assert(order == 0); % only implemented now

    function [py pv] = sparsePrior(x)
        x = reshape(x,(CSP-dim),L);
        
        py = 0;
        pv = zeros(CSP,L);
        
        for l = 1:L
            for si = 1:R
                alsi = x(dim*(si-1)+(1:dim),l);
                
%                 lambdaspsi = lambdasp;
                lambdaspsi = lambdasp*(max(scales)^alpha/(scales(si)^alpha));
                
                if norm(alsi) > c && lambdaspsi > 0.0
                    py = py + lambdaspsi*(max(log(norm(alsi)),log(c))-log(c));
                    pv(dim*(si-1)+dim+(1:dim),l) = pv(dim*(si-1)+dim+(1:dim),l)+lambdaspsi*alsi/norm(alsi)^2;
%                     py = py + lambdaspsi*norm(alsi);
%                     pv(dim*(si-1)+dim+(1:dim),l) = pv(dim*(si-1)+dim+(1:dim),l)+lambdaspsi*alsi/norm(alsi);
                end
            end
        end
        
        pv = reshape(pv,CSP*L,1);
    end

prior = @sparsePrior;

end