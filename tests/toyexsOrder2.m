% produce toy figures for hok paper

% LDDMM options
clear lddmmoptions
lddmmoptions.scale = 1.0; % Gaussian kernels
lddmmoptions.energyweight = [1 12]; % weighting between energy of curve and match
lddmmoptions.energyweight = lddmmoptions.energyweight/sum(lddmmoptions.energyweight);
lddmmoptions.order = 2;

% output options
% global globalOptions;
% globalOptions.verbose = true;
visoptions.dim = 2;
visoptions.Ngridpoints = 25;
visoptions.margin = 1.5;
visoptions.skip = 4;
visoptions.fixedDim = true;
visoptions.minx = 3.5;
visoptions.maxx = 6.5;
visoptions.miny = 3.5;
visoptions.maxy = 6.5;
% scale relative to grid:
scaleInGridUnits = lddmmoptions.scale/(2*visoptions.margin/(visoptions.Ngridpoints-1))

options = getDefaultOptions();

% one point examples
dim = 2;
moving = [
    5.0 5.0;
]';    
% fixed - matching against these
fixed = [
     5.0 5.0 reshape(eye(dim),dim^2,1)';
]';
[methods lddmmoptions1] = setupPointLDDMM(moving,fixed,[],lddmmoptions);
visualizer = getPointVisualizer(methods.transport,moving,[],visoptions);
% expansion
figure(1)
A1 = zeros(dim,dim);
% A1 = eye(dim);
A2 = zeros(dim,dim,dim);
% A2(1,2,1) = 1; A2(1,1,2) = 1;
% A2(1,1,1) = 1;
A2(1,2,2) = 1;
A2 = 0.2*A2;
x = [0 0 reshape(A1,dim^2,1)' reshape(A2,dim^3,1)']';
visualizer(x);
[res,Gt] = methods.transport(x,moving);
E = methods.pathEnergy(x,Gt);
E
Gt.x, diff(Gt.x)

% two point examples
dim = 2;
moving = [
    -1.0 0.0;
    1.0 0.0;
]';    
[methods lddmmoptions1] = setupPointLDDMM(moving,[],[],lddmmoptions);
visualizer = getPointVisualizer(methods.transport,moving,[],visoptions);
% expansion
figure(1)
A1 = zeros(dim,dim);
% A1 = eye(dim);
A2 = zeros(dim,dim,dim);
% A2(1,2,1) = 1; A2(1,1,2) = 1;
A2(1,1,1) = 1;
% A2(1,2,2) = 1;
x = [0 0 reshape(A1,dim^2,1)' reshape(A2,dim^3,1)'; 
     0 0 reshape(A1,dim^2,1)' reshape(A2,dim^3,1)']';
visualizer(x);
[res,Gt] = methods.transport(x,moving);
E = methods.pathEnergy(x,Gt);
E
Gt.x, diff(Gt.x)
