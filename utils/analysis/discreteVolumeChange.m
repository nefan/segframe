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

function vChange = discreteVolumeChange(x,IF,transport)

% obtain moving image result
[grid1{1} grid1{2} grid1{3}] = meshgrid(1:size(IF,2),1:size(IF,1),0);

% transport
grid0 = transport(x,grid1,true);

vChange = 0;
for i = 1:size(IF,1)-1
    for j = 1:size(IF,2)-1
        vx = [grid0{1}(j,i+1)-grid0{1}(j,i); grid0{2}(j,i+1)-grid0{2}(j,i)];
        vy = [grid0{1}(j+1,i)-grid0{1}(j,i); grid0{2}(j+1,i)-grid0{2}(j,i)];
        
        vChange = vChange + abs((abs(det([vx vy]))-1));
    end
end
vChange = vChange/((size(IF,1)-1)*(size(IF,2)-1));

end