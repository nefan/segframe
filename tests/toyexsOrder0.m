% produce toy figures (originally for higher-order momenta paper)

% LDDMM options
clear lddmmoptions
lddmmoptions.scale = 1.0; % Gaussian kernels
lddmmoptions.energyweight = [1 12]; % weighting between energy of curve and match
lddmmoptions.energyweight = lddmmoptions.energyweight/sum(lddmmoptions.energyweight);
lddmmoptions.order = 0;

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
     5.0 5.0;
]';
[methods lddmmoptions1] = setupPointLDDMM(moving,fixed,[],lddmmoptions);
visualizer = getPointVisualizer(methods.transport,moving,[],visoptions);
% translation
figure(1)
x = [1.25 0]';
visualizer(x);
% approximated derivatives
figure(8)
moving = [
    5.25 5.0;
    4.75 5.0;
    5.0 5.25;
    5.0 4.75;
]';    
% fixed
fixed = [
     5.2 5.0;
     4.8 5.0;
     5.0 5.2;
     5.0 4.8;
     ]';
[methods lddmmoptions1] = setupPointLDDMM(moving,fixed,[],lddmmoptions);
visualizer = getPointVisualizer(methods.transport,moving,[],visoptions);
x = [
     2.3 0;
     -2.3 0;
     0 2.3;
     0 -2.3;
     ]';
visualizer(x);


% two points examples
dim = 2;
moving = [
    4.5 5.0;
    5.5 5.0;
]';    
% fixed - matching against these
fixed = [
     4.0 5.0;
     5.0 5.0;
]';
[methods lddmmoptions1] = setupPointLDDMM(moving,fixed,[],lddmmoptions);
visualizer = getPointVisualizer(methods.transport,moving,[],visoptions);
% rotation
figure(4)
x = 10*[
    1 0;
    -1 0;
    ]';
x = reshape(x,[],1);
tend = 1.0;
visualizer(x,tend);
[res,Gt] = methods.transport(x,moving,false,tend);
E = methods.pathEnergy(x,Gt,tend);
% E
% Gt.x, diff(Gt.x)
