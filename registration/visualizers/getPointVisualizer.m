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

function visualizer = getPointVisualizer(transport,moving,fixed,visoptions)
% intial points: red
% end points: green
% results: black

dim = visoptions.dim;
% parameters
Ngridpoints = getOption(visoptions,'Ngridpoints',16);
margin = getOption(visoptions,'margin',1);
skip = getOption(visoptions,'skip',1); % only display grid for every skip gridpoints
fixedDim = getOption(visoptions,'fixedDim',false);
if fixedDim
    mminx = getOption(visoptions,'minx',-1);
    mmaxx = getOption(visoptions,'maxx',-1);
    mminy = getOption(visoptions,'miny',-1);
    mmaxy = getOption(visoptions,'maxy',-1);
end


    function pointVisualizer(x,varargin)
        tend = 1;
        if size(varargin,2) > 0
            tend = varargin{1};
        end

        
        % determine area
        res = transport(x,moving,false,tend);
        if ~isempty(fixed)
            minx = min([moving(1,:) fixed(1,:) res(1,:)]);
            maxx = max([moving(1,:) fixed(1,:) res(1,:)]);
            miny = min([moving(2,:) fixed(2,:) res(2,:)]);
            maxy = max([moving(2,:) fixed(2,:) res(2,:)]);        
        else
            minx = min([moving(1,:) res(1,:)]);
            maxx = max([moving(1,:)]);
            miny = min([moving(2,:) res(2,:)]);
            maxy = max([moving(2,:) res(2,:)]);        
        end
        if dim == 3
            assert(false); % not implemented
        end
        
        % create grid
        if ~fixedDim
            mminx = minx-margin; mmaxx = maxx+margin;
            mminy = miny-margin; mmaxy = maxy+margin;
        end
        if dim == 3
            assert(false);
        end
        
        if dim == 2
            [grid0{1} grid0{2} grid0{3}] = meshgrid(...
                mminx:(mmaxx-mminx)/((Ngridpoints-1)*skip):mmaxx,...
                mminy:(mmaxy-mminy)/((Ngridpoints-1)*skip):mmaxy,...
                0);   
        else
            [grid0{1} grid0{2} grid0{3}] = meshgrid(...
                mminx:(mmaxx-mminx)/((Ngridpoints-1)*skip):mmaxx,...
                mminy:(mmaxy-mminy)/((Ngridpoints-1)*skip):mmaxy,...
                mminz:(mmaxz-mminz)/((Ngridpoints-1)*skip):mmaxz);               
        end
        
        printStatus('integrating grid');
        pgrid0 = [reshape(grid0{1},1,[]); reshape(grid0{2},1,[]); reshape(grid0{3},1,[])];
        pgrid1 = transport(x,pgrid0,false,tend);
        grid1{1} = reshape(pgrid1(1,:),size(grid0{1})); 
        grid1{2} = reshape(pgrid1(2,:),size(grid0{1}));  
        grid1{3} = zeros(size(grid0{1}));
        
        % show grid
        clf
        show2dgrid(grid1,'r',skip,grid0,visoptions);        
        showpoints(moving,'r',visoptions);
        showpoints(res,'g',visoptions);
        if ~isempty(fixed)
            showpoints(fixed(1:dim,:),'k',visoptions);
        end
    end

visualizer = @pointVisualizer;

end