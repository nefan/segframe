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

function show2dgrid(grid,color,skip,grid0,visoptions)

dim = visoptions.dim;
doColor = getOption(visoptions,'colorGrid');

axis off
axis equal

lgrid1 = squeeze(grid{1});
lgrid2 = squeeze(grid{2});
lgrid3 = squeeze(grid{3});

if ~doColor
    mesh(lgrid1(1:skip:end,1:skip:end),lgrid2(1:skip:end,1:skip:end),lgrid3(1:skip:end,1:skip:end),'EdgeColor','k','FaceAlpha',0.0,'LineWidth',1.0);
else    
    % 2d quick hack for coloring with straing tensor
    % done shortly before IPMI deadline :-)

    % strain tensor colors
    lgrid01 = squeeze(grid0{1});
    lgrid02 = squeeze(grid0{2});
    lgrid03 = squeeze(grid0{3});
    gridcolors = zeros(size(lgrid1));
    for i = 1:size(lgrid1,1) % vertical
        for j = 1:size(lgrid1,2) % horizontal
            p = [lgrid1(i,j) lgrid2(i,j)];

            xbar = [lgrid01(i,j) lgrid02(i,j)]';
            ybar = p';
            xi = [];
            yi = [];
            if j > 1
                xi = [xi [lgrid01(i,j-1) lgrid02(i,j-1)]'];
                yi = [yi [lgrid1(i,j-1) lgrid2(i,j-1)]'];
            end
            if j < size(lgrid1,2)
                xi = [xi [lgrid01(i,j+1) lgrid02(i,j+1)]'];
                yi = [yi [lgrid1(i,j+1) lgrid2(i,j+1)]'];
            end
            if i > 1
                xi = [xi [lgrid01(i-1,j) lgrid02(i-1,j)]'];
                yi = [yi [lgrid1(i-1,j) lgrid2(i-1,j)]'];
            end
            if i < size(lgrid1,2)
                xi = [xi [lgrid01(i+1,j) lgrid02(i+1,j)]'];
                yi = [yi [lgrid1(i+1,j) lgrid2(i+1,j)]'];
            end
            sigmaxy = zeros(2,2);
            sigmaxx = zeros(2,2);
            for k = 1:size(xi,2)
                sigmaxy = sigmaxy+1/size(xi,2)*(xi(:,k)-xbar)*(yi(:,k)-ybar)';
                sigmaxx = sigmaxx+1/size(xi,2)*(xi(:,k)-xbar)*(xi(:,k)-xbar)';
            end
            A = sigmaxy'*inv(sigmaxx);
            gridcolors(i,j) = log(trace(A*A'));
        end
    end
    % remember maximal values for several runs of scalematch
%     global gridmin;
%     global gridmax;
    gridmin = []; gridmax = [];
    if gridmin
        gridmin = min(min(min(gridcolors)),gridmin);
    else
        gridmin = min(min(gridcolors));
    end
    if gridmax
        gridmax = max(max(max(gridcolors)),gridmax);
    else
        gridmax = max(max(gridcolors));
    end
    gridcolors = (gridcolors-gridmin)/(gridmax-gridmin);

    for i = 1:size(lgrid1,1) % vertical
        for j = 1:size(lgrid1,2) % horizontal
            p = [lgrid1(i,j) lgrid2(i,j)];        
            colorv = max(gridcolors(i,j),0);

            cm = colormap();
            color = cm(floor(colorv*(size(colormap,1)-1))+1,:);

            if j > 1 && mod(i-1,skip) == 0
                line([p(1) lgrid1(i,j-1)],[p(2) lgrid2(i,j-1)],'Color',color,'LineWidth',1.2);
            end
            if i > 1 && mod(j-1,skip) == 0
                line([p(1) lgrid1(i-1,j)],[p(2) lgrid2(i-1,j)],'Color',color,'LineWidth',1.2);
            end
        end
    end

end