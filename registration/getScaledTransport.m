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

function transport = getScaledTransport(unscaledTransport,scale)

    function res = ltransport(x, points, varargin)
        
        if iscell(points)
            g0 = points;
            g0{1} = scale*g0{1};
            g0{2} = scale*g0{2};
            g0{3} = scale*g0{3};       
        else
            g0{1} = scale*points(1,:);
            g0{2} = scale*points(2,:);
            if dim == 2
                g0{3} = zeros(size(g0{1}));
            else
                g0{3} = scale*points(3,:);
            end
        end
        
        g1 = unscaledTransport(x, g0, varargin{:});
        
        if iscell(points)
            g1{1} = 1/scale*g1{1};
            g1{2} = 1/scale*g1{2};
            g1{3} = 1/scale*g1{3};
            res = g1;
        else
            if dim == 2
                res = 1/scale*[g1{1}; g1{2}];
            else
                res = 1/scale*[g1{1}; g1{2}; g1{3}];
            end
        end
    end

transport = @ltransport;

end