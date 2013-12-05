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
imageoptions.scale = getOption(imageoptions,'scale',lddmmoptions.scale);
imageoptions.range = ceil(getOption(imageoptions,'range',3*imageoptions.scale));
imageoptions.margin = ceil(getOption(imageoptions,'margin',0));
imageoptions.background = getOption(imageoptions,'background',0);
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
    spacingX = (size(IF,1)-2*imageoptions.margin)/lddmmoptions.NpointsX;
    spacingY = (size(IF,2)-2*imageoptions.margin)/lddmmoptions.NpointsY;

    [grid0{1} grid0{2}] = meshgrid(spacingX/2+imageoptions.margin:spacingX:size(IF,1)-imageoptions.margin,spacingY/2+imageoptions.margin:spacingY:size(IF,2)-imageoptions.margin);

    lddmmoptions.L = lddmmoptions.NpointsX*lddmmoptions.NpointsY;
    moving = [reshape(grid0{1},1,lddmmoptions.L); reshape(grid0{2},1,lddmmoptions.L)];
    imageoptions.h = max(spacingX,spacingY);
end
lddmmoptions.L = size(moving,2);
imageoptions.moving = moving;


% lddmm stuff
[methods lddmmoptions optimoptions] = setupLDDMM(moving,lddmmoptions,varargin{:});

% options from lddmm
imageoptions.dimq = lddmmoptions.dimq;

% methods
switch imageoptions.measure
    case 'LOI'
        imageU = getImageULOI(IM,IF,moving,imageoptions,lddmmoptions.dim,lddmmoptions.L);
    case 'L2'
        imageU = getImageUL2(IM,IF,moving,imageoptions,lddmmoptions.dim,lddmmoptions.L);
end        
        
methods.F = getPointLDDMMF(imageU,moving,...
    methods,lddmmoptions);

end