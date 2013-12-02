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

function [methods lddmmoptions optimoptions] = setupLDDMM(moving,lddmmoptions,varargin)

optimoptions = [];
if size(varargin,2) > 0
    optimoptions = varargin{1};
end

% dim
assert(lddmmoptions.dim == 2 || lddmmoptions.dim == 3);
lddmmoptions.cdim = 3; % all computations done in dim 3

% scales
assert(isfield(lddmmoptions,'scale') || isfield(lddmmoptions,'scales'));
if isfield(lddmmoptions,'scale')
    assert(~isfield(lddmmoptions,'scales') || length(lddmmoptions.scales) == 1);
    lddmmoptions.R = 1;
    lddmmoptions.scales = [lddmmoptions.scale];
    lddmmoptions.scaleweight = [1];
else
    lddmmoptions.R = length(lddmmoptions.scales);
    lddmmoptions.scaleweight = ones(lddmmoptions.R);
end

% order (> 0 implies higher order kernels)
if ~isfield(lddmmoptions,'order')
    lddmmoptions.order = 0;
end
assert(lddmmoptions.order >= 0 || lddmmoptions.order <= 2);
assert(lddmmoptions.order ==0 || lddmmoptions.R == 1); % scale currently only supported for order 0

% data structures
switch lddmmoptions.order
    case 0
        lddmmoptions.dimq = lddmmoptions.dim;
        lddmmoptions.CSP = lddmmoptions.dim+lddmmoptions.dim*lddmmoptions.R; % column size particles
        lddmmoptions.cCSP = lddmmoptions.cdim+lddmmoptions.cdim*lddmmoptions.R; % for computations
    case 1
        lddmmoptions.dimq = lddmmoptions.dim+lddmmoptions.dim^2;
        lddmmoptions.CSP = lddmmoptions.dim*(2+2*lddmmoptions.dim); % forward column size particles
        lddmmoptions.cCSP = lddmmoptions.cdim*(2+2*lddmmoptions.cdim); % for computations
    case 2
        lddmmoptions.dimq = lddmmoptions.dim+lddmmoptions.dim^2+lddmmoptions.dim^3;
        lddmmoptions.CSP = 2*lddmmoptions.dim+2*lddmmoptions.dim^2+2*lddmmoptions.dim^3; % forward column size particles
        lddmmoptions.cCSP = 2*lddmmoptions.cdim+2*lddmmoptions.cdim^2+2*lddmmoptions.cdim^3; % for computations
end

% tolerances
% optimoptions.tolFun = 1e-7;
% optimoptions.tolX = 1e-7;
lddmmoptions.epsilon = 10e-10;

% misc
if ~isfield(optimoptions,'verbose')
    optimoptions.verbose = true;
end
% debugging
% optimoptions.derivativeCheck = true;
% lddmmoptions.testC = true;

% methods
if getOption(lddmmoptions,'sparsity')
    methods.prior = getSparsePrior(lddmmoptions);
end
methods.optimizer = getDefaultOptimizer(optimoptions);
methods.getInitialData = getPointLDDMMInitialData(moving,lddmmoptions);
methods.pointPath = getPointPath(lddmmoptions);
methods.gradTransport = getPointGradTransport(lddmmoptions);
methods.transport = getPointLDDMMTransport(getPointPath(lddmmoptions),moving,lddmmoptions);
methods.pathEnergy = getPointPathEnergy(lddmmoptions);
methods.getMoving = @() moving;

end