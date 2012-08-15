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

function visualizer = getPointLDDMMVisualizer(transport,moving,fixed,visoptions,lddmmoptions)
% intial points: red
% end points: black
% match points: green

dim = lddmmoptions.dim;
pointVisualizer = getPointVisualizer(transport,moving,fixed,visoptions);

    function pointLDDMMVisualizer(x)
        
        pointVisualizer(x);
        showmomentum(reshape([moving; reshape(x,[],lddmmoptions.L)],[],lddmmoptions.L),visoptions,lddmmoptions);

    end

visualizer = @pointLDDMMVisualizer;

end