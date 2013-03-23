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
dim = lddmmoptions.dim;
cdim = lddmmoptions.cdim;

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
assert(lddmmoptions.order >= 0 || lddmmoptions.order <= 1);
if lddmmoptions.order == 0
    lddmmoptions.dimX1 = (1+lddmmoptions.R)*dim; % x dimension
else
    lddmmoptions.dimX1 = (1+dim*lddmmoptions.R)*dim;
end

% data structures
if lddmmoptions.order == 0
    lddmmoptions.CSP = dim*(1+lddmmoptions.R); % column size particles
    lddmmoptions.cCSP = cdim*(1+lddmmoptions.R); % for computations
    lddmmoptions.iCSP = lddmmoptions.R*dim; % input
    
    % indexing
    lddmmoptions.Iphi = @() 1:cdim;
    lddmmoptions.Irho = @(s) cdim*s+1:cdim;
    
else
    lddmmoptions.CSP = dim+lddmmoptions.R*(dim+dim^2); % obsolete
    lddmmoptions.iCSP = lddmmoptions.R*(dim+dim^2); % input
    lddmmoptions.fCSP = cdim+cdim^2+lddmmoptions.R*cdim; % forward
    lddmmoptions.bCSP = cdim+cdim^2+lddmmoptions.R*(cdim+cdim^2); % backwards
    
    % indexing
    phioffset = 0;
    Dphioffset = phioffset+cdim;
    muoffset = Dphioffset+cdim*cdim;
    mujoffset = muoffset+cdim;
    iscaleoffset = dim+dim^2;
    fscaleoffset = cdim;
    bscaleoffset = cdim+cdim^2;

    % input
    lddmmoptions.iImu = @(s) (s-1)*iscaleoffset+(1:dim);    
    lddmmoptions.iImuj = @(s) dim+(s-1)*iscaleoffset+(1:dim^2);    
    
    % forward
    lddmmoptions.fIphi = @() phioffset+(1:cdim);
    lddmmoptions.fIDphi = @() Dphioffset+(1:cdim^2);
    lddmmoptions.fIDphic = @(c) Dphioffset+cdim*(c-1)+(1:cdim);
    lddmmoptions.fImu = @(s) muoffset+(s-1)*fscaleoffset+(1:cdim);

    % derivatives
    lddmmoptions.bIphi = @() phioffset+(1:cdim);
    lddmmoptions.bIDphi = @() Dphioffset+(1:cdim^2);
    lddmmoptions.bImu = @(s) muoffset+(s-1)*bscaleoffset+(1:cdim);    
    lddmmoptions.bImu = @(s) mujoffset+(s-1)*bscaleoffset+(1:cdim^2);
end

% tolerances
% optimoptions.tolFun = 1e-7;
% optimoptions.tolX = 1e-7;
lddmmoptions.epsilon = 10e-10;

% misc
if ~isfield(optimoptions,'verbose')
    optimoptions.verbose = true;
end
% optimoptions.derivativeCheck = true;
% lddmmoptions.testC = true; % for dim 3

% methods
if isfield(lddmmoptions,'reverse')
    assert(~lddmmoptions.reverse); % not implemented
end
if getOption(lddmmoptions,'sparsity')
    methods.prior = getSparsePrior(lddmmoptions);
end
methods.optimizer = getDefaultOptimizer(optimoptions);
methods.getInitialData = getPointLDDMMInitialData(moving,lddmmoptions);
methods.pointPath = getPointPath(moving,lddmmoptions);
methods.gradTransport = getPointGradTransport(lddmmoptions);
methods.transport = getPointLDDMMTransport(getPointPath(moving,lddmmoptions),lddmmoptions);
methods.pathEnergy = getPointPathEnergy(lddmmoptions);

end