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
imageoptions.measure = 'L2';

% LDDMM options
clear lddmmoptions
lddmmoptions.energyweight = [1 2^4]; % weighting between energy of curve and match
lddmmoptions.energyweight = lddmmoptions.energyweight/sum(lddmmoptions.energyweight);
lddmmoptions.order = 0;

% control points
lddmmoptions.NpointsX = 2;
lddmmoptions.NpointsY = 2;

% output options
% global globalOptions;
% globalOptions.verbose = true;
clear visoptions
visoptions = [];

% data
DD = load('tests/data/images1.mat');

match = [];
IM = DD.Im;
IF = DD.If;
% scale
spacingX = size(IF,1)/lddmmoptions.NpointsX;
spacingY = size(IF,2)/lddmmoptions.NpointsY;
imageoptions.scale = max(spacingX,spacingY)/2;
lddmmoptions.scale = imageoptions.scale; % Gaussian kernels

options = getDefaultOptions();

% optimization
optimoptions.numDiff = true;

[methods lddmmoptions imageoptions] = setupImageLDDMM(IM,IF,imageoptions,lddmmoptions,optimoptions);

% visualizer
visualizer = getImageVisualizerL2(methods.transport,IM,IF,visoptions,imageoptions);
methods.iterationVisualizer = visualizer;

result = runRegister(methods, options);

% visualize
visualizer(result);

figure(11), clf
visoptions.dim = lddmmoptions.dim;
visoptions.margin = 0;
pointVisualizer = getPointVisualizer(methods.transport,methods.getMoving(),[],visoptions);
pointVisualizer(result);