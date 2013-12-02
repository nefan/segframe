% produce toy figures for hok paper

% LDDMM options
clear lddmmoptions
lddmmoptions.scale = 1.0; % Gaussian kernels
lddmmoptions.energyweight = [1 12]; % weighting between energy of curve and match
lddmmoptions.energyweight = lddmmoptions.energyweight/sum(lddmmoptions.energyweight);
lddmmoptions.order = 1;

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
x = [0 0 reshape(eye(dim),dim^2,1)']';
visualizer(x);
% contraction
figure(2)
x = [0 0 reshape(-eye(dim),dim^2,1)']';
visualizer(x);
% rotation
figure(3)
v = -4/8*pi;
x = [0 0 reshape([cos(v) sin(v); -sin(v) cos(v)],dim^2,1)']';
visualizer(x);
% translations
figure(7)
x = [1.25 0 zeros(dim^2,1)']';
visualizer(x);
% approximated derivatives
figure(8)
moving = [
    5.25 5.0;
    4.75 5.0;
    5.0 5.25;
    5.0 4.75;
]';    
% fixed - matching against these
fixed = [
     5.2 5.0 reshape(eye(dim),dim^2,1)';
     4.8 5.0 reshape(eye(dim),dim^2,1)';
     5.0 5.2 reshape(eye(dim),dim^2,1)';
     5.0 4.8 reshape(eye(dim),dim^2,1)';
     ]';
[methods lddmmoptions1] = setupPointLDDMM(moving,fixed,[],lddmmoptions);
visualizer = getPointVisualizer(methods.transport,moving,[],visoptions);
x = [
     2.3 0 zeros(dim^2,1)';
     -2.3 0 zeros(dim^2,1)';
     0 2.3 zeros(dim^2,1)';
     0 -2.3 zeros(dim^2,1)';]';
visualizer(x);


% two points examples
dim = 2;
moving = [
    3.5 5.0;
    6.5 5.0;
]';    
% fixed - matching against these
fixed = [
     4.0 5.0 reshape(eye(dim),dim^2,1)';
     6.0 5.0 reshape(eye(dim),dim^2,1)';
]';
[methods lddmmoptions1] = setupPointLDDMM(moving,fixed,[],lddmmoptions);
visualizer = getPointVisualizer(methods.transport,moving,[],visoptions);
% rotation
figure(4)
v = 4/8*pi;
x = [
    0 0 reshape([cos(v) sin(v); -sin(v) cos(v)],dim^2,1)'
    0 0 reshape([cos(v) sin(v); -sin(v) cos(v)],dim^2,1)'
    ]';
visualizer(reshape(x,[],1));

% two point match examples
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
[methods lddmmoptions1] = setupPointLDDMM(moving,fixed,[],lddmmoptions);
visualizer = getPointVisualizer(methods.transport,moving,fixed,visoptions);
figure(5)
result = runRegister(methods, options);
visualizer(result);
% two point match examples
dim = 2;
moving = [
    3.5 4.0;
    6.5 4.0;
]';    
% fixed - matching against these
v = -4/8*pi;
fixed = [
     4.0 6.0 reshape([cos(v) sin(v); -sin(v) cos(v)],dim^2,1)';
     6.0 6.0 reshape([cos(-v) sin(-v); -sin(-v) cos(-v)],dim^2,1)';
]';
[methods lddmmoptions1] = setupPointLDDMM(moving,fixed,[],lddmmoptions);
visualizer = getPointVisualizer(methods.transport,moving,fixed,visoptions);
figure(6)
result = runRegister(methods, options);
visualizer(result);
