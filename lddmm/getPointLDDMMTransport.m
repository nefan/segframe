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
scales = lddmmoptions.scales;
scaleweight = lddmmoptions.scaleweight;
energyweight = lddmmoptions.energyweight;

    function res = ltransport(x, points, varargin)
        backwards = false;
        if size(varargin,2) > 0
            backwards = varargin{1};
        end
        
        ic = pathIC(x,moving,lddmmoptions);        
        Gt = pointPath(ic); % rhot is a solver structure
        
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
        if order == 0
            g1 = shootgrid(g0,Gt,lddmmoptions,backwards);
        else
            g1 = shootgridDKernels(x,g0,Gt,lddmmoptions,backwards);
        end
        if iscell(points)
            res = g1;
        else
            if dim == 2
                res = [g1{1}; g1{2}];
            else
                res = [g1{1}; g1{2}; g1{3}];
            end
        end
    end

transport = @ltransport;

end