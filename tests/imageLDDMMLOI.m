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

% control points
lddmmoptions.NpointsX = 10;
lddmmoptions.NpointsY = 10;

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

% % pre-spline
% [grid{1} grid{2} grid{3}] = meshgrid(1:size(IF,1),1:size(IF,2),0);
% IM = reshape(Spline2D([reshape(grid{1},numel(IF),1) reshape(grid{2},numel(IF),1)],...
%     zeros(numel(IF),1),IM,[0 0],double([1 1]),1),size(IM))';
% IF = reshape(Spline2D([reshape(grid{1},numel(IF),1) reshape(grid{2},numel(IF),1)],...
%     zeros(numel(IF),1),IF,[0 0],double([1 1]),1),size(IF))';

options = getDefaultOptions();

[methods lddmmoptions] = setupImageLDDMM(IM,IF,imageoptions,lddmmoptions);

% visualizer
visualizer = getImageVisualizer(methods.transport,IM,IF,visoptions,imageoptions);
methods.iterationVisualizer = visualizer;

result = runRegister(methods, options);

% visualize
visualizer(result);

figure(11), clf
visoptions.dim = lddmmoptions.dim;
visoptions.margin = 0;
pointVisualizer = getPointVisualizer(methods.transport,methods.getMoving(),[],visoptions);
pointVisualizer(result);