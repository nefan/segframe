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

function [methods lddmmoptions imageoptions] = setupImageLDDMM(IM,IF,imageoptions,lddmmoptions,varargin)

% image specific
imageoptions.reverse = false; % if true, measure in moving image
lddmmoptions.reverse = ~imageoptions.reverse;
imageoptions.scale = lddmmoptions.scale;
imageoptions.range = getOption(imageoptions,'range',3*imageoptions.scale);
% order: 0: only points information, 1: point and derivative information
if ~isfield(lddmmoptions,'order')
    lddmmoptions.order = 0;
end
imageoptions.order = lddmmoptions.order;

% dim
imageoptions.dim = length(size(IF));
lddmmoptions.dim = imageoptions.dim;

% points
isequal(size(IM),size(IF));
if isfield(lddmmoptions,'controlPoints')
    moving = lddmmoptions.controlPoints;
else
    assert(imageoptions.dim == 2);
    spacingX = size(IF,1)/lddmmoptions.NpointsX;
    spacingY = size(IF,2)/lddmmoptions.NpointsY;

    [grid0{1} grid0{2}] = meshgrid(spacingX/2:spacingX:size(IF,1),spacingY/2:spacingY:size(IF,2));

    lddmmoptions.L = lddmmoptions.NpointsX*lddmmoptions.NpointsY;
    moving = [reshape(grid0{1},1,lddmmoptions.L); reshape(grid0{2},1,lddmmoptions.L)]; 
end
lddmmoptions.L = size(moving,2);


% lddmm stuff
[methods lddmmoptions optimoptions] = setupLDDMM(moving,lddmmoptions,varargin{:});

% methods
methods.F = getPointLDDMMF(getImageU(IM,IF,moving,imageoptions,lddmmoptions.dim,lddmmoptions.L),...
    methods,lddmmoptions);


end