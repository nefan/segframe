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

function M = mmap(f,varargin)
%
% creating matrix M with dimensions dim by mapping function f
% on all index tuples
%

if size(varargin,2) == 1
    [d1]=ndgrid(1:varargin{1});
    M = arrayfun(f,d1);
else if size(varargin,2) == 2
    [d1,d2]=ndgrid(1:varargin{1},1:varargin{2});
    M = arrayfun(f,d1,d2);
else if size(varargin,2) == 3
    [d1,d2,d3]=ndgrid(1:varargin{1},1:varargin{2},1:varargin{3});
    M = arrayfun(f,d1,d2,d3);
else if size(varargin,2) == 4
    [d1,d2,d3,d4]=ndgrid(1:varargin{1},1:varargin{2},1:varargin{3},1:varargin{4});
    M = arrayfun(f,d1,d2,d3,d4);
else if size(varargin,2) == 5
    [d1,d2,d3,d4,d5]=ndgrid(1:varargin{1},1:varargin{2},1:varargin{3},1:varargin{4},1:varargin{5});
    M = arrayfun(f,d1,d2,d3,d4,d5);
else if size(varargin,2) == 6
    [d1,d2,d3,d4,d5,d6]=ndgrid(1:varargin{1},1:varargin{2},1:varargin{3},1:varargin{4},1:varargin{5},1:varargin{6});
    M = arrayfun(f,d1,d2,d3,d4,d5,d6);
else if size(varargin,2) == 7
    [d1,d2,d3,d4,d5,d6,d7]=ndgrid(1:varargin{1},1:varargin{2},1:varargin{3},1:varargin{4},1:varargin{5},1:varargin{6},1:varargin{7});
    M = arrayfun(f,d1,d2,d3,d4,d5,d6,d7);    
else if size(varargin,2) == 8
    [d1,d2,d3,d4,d5,d6,d7,d8]=ndgrid(1:varargin{1},1:varargin{2},1:varargin{3},1:varargin{4},1:varargin{5},1:varargin{6},1:varargin{7},1:varargin{8});
    M = arrayfun(f,d1,d2,d3,d4,d5,d6,d7,d8);    
else if size(varargin,2) > 8
        assert(false);
else
    M = [];
end
end
end
end
end
end
end
end
end
end