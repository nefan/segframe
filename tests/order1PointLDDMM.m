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

% test for order 1point LDDMM registration

% LDDMM options
clear lddmmoptions
lddmmoptions.scale = 1.5; % Gaussian kernels
lddmmoptions.energyweight = [1 8]; % weighting between energy of curve and match
lddmmoptions.energyweight = lddmmoptions.energyweight/sum(lddmmoptions.energyweight);
lddmmoptions.order = 1;

% output options
% global globalOptions;
% globalOptions.verbose = true;
clear visoptions
visoptions.dim = 2;

% two point match example
dim = 2;
moving = [
    3.5 4.0;
    6.5 4.0;
]';    
% fixed - matching against these
fixed = [
     4.0 6.0 reshape(2*eye(dim),dim^2,1)';
     6.0 6.0 reshape(2*eye(dim),dim^2,1)';
]';

options = getDefaultOptions();

[methods lddmmoptions1] = setupPointLDDMM(moving,fixed,[],lddmmoptions);

result = runRegister(methods, options);

% visualize
visualizer = getPointVisualizer(methods.transport,moving,fixed,visoptions);
visualizer(result);
