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

% test for LDDMM image registration

% image options
% imageoptions.samplesPerPoint = 100;
imageoptions.randomSampling = false;

% LDDMM options
clear lddmmoptions
lddmmoptions.scale = 25; % Gaussian kernels
lddmmoptions.energyweight = [1 2^8]; % weighting between energy of curve and match
lddmmoptions.energyweight = lddmmoptions.energyweight/sum(lddmmoptions.energyweight);
lddmmoptions.order = 1;

% control points
% lddmmoptions.controlPoints = [130; 130];
lddmmoptions.NpointsX = 10;
lddmmoptions.NpointsY = 10;

% output options
% global globalOptions;
% globalOptions.verbose = true;
clear visoptions
visoptions.dim = 2;

% data
DD = load('tests/data/images1.mat');

match = [];
IM = DD.Im;
IF = DD.If;

options = getDefaultOptions();

[methods lddmmoptions imageoptions] = setupImageLDDMM(IM,IF,imageoptions,lddmmoptions);

result = runRegister(methods, options);

% visualize
visualizer = getImageVisualizer(methods.transport,IM,IF,visoptions,imageoptions);
visualizer(result);

figure(11), clf
visoptions.dim = lddmmoptions.dim;
visoptions.margin = 0;
pointVisualizer = getPointVisualizer(methods.transport,methods.getMoving(),[],visoptions);
pointVisualizer(result);