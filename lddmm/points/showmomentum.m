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

function showmomentum(x,visoptions,lddmmoptions)

dim = lddmmoptions.dim;
L = lddmmoptions.L;
R = lddmmoptions.R;

dispScale = getOption(visoptions,'dispScale');
if isempty(dispScale)
    dispScale = 5;
end

rho = reshape(x,(1+R)*dim,L);

hold on

for l = 1:L
    % point
    p = rho(:,l);
    
    % value
    for s = 1:R
        color = [0 s/3 s/3];
        if dim == 2
            if norm(p(dim*(s-1)+dim+(1:dim))) < 0.035 % stuff for sparsity plots
            else
            if norm(p(dim*(s-1)+dim+(1:dim))) < 0.5
                quiver(p(1),p(2),p(dim*(s-1)+1+dim),p(dim*(s-1)+2*dim),'Color',color,'LineWidth',3,'MaxHeadSize',1.3);
            else
                quiver(p(1),p(2),p(dim*(s-1)+1+dim),p(dim*(s-1)+2*dim),'Color',color,'LineWidth',3,'MaxHeadSize',0.4);
            end
            end
        else
            quiver3(p(1),p(2),p(3),p(dim*(s-1)+1+dim),p(dim*(s-1)+1+dim+1),p(dim*(s-1)+2*dim),'Color',color,'LineWidth',2);
        end
    end
    
%     % scale
%     circle = [];
%     r = p(5);
%     for t = 0:.01:2*pi
%         circle = [circle [p(1)+r*cos(t) p(2)+r*sin(t)]'];
%     end
%     plot(circle(1,:),circle(2,:),'.','MarkerSize',1,'MarkerEdgeColor','r');    
end
