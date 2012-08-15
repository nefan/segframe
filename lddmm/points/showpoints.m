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

function showpoints(x,color,visoptions)

dim = visoptions.dim;

dispScale = getOption(visoptions,'dispScale');
if isempty(dispScale)
    dispScale = 5;
end

x = reshape(x,dim,[]);
x = x(1:dim,:);

hold on

for l = 1:size(x,2)
    % point
    p = x(:,l);
    if dim == 2
        plot(p(1),p(2),'.','MarkerSize',dispScale*4.2,'MarkerFaceColor',color,'MarkerEdgeColor',color);
    else
        plot3(p(1),p(2),p(3),'.','MarkerSize',dispScale*4.2,'MarkerFaceColor',color,'MarkerEdgeColor',color);
    end
end