% register binary images

% random numbers
% s = RandStream.create('mt19937ar','seed',342484);
% s = RandStream.create('mt19937ar','seed',3853873);
% RandStream.setDefaultStream(s);

% LDDMM options
clear lddmmoptions
lddmmoptions.scale = 50; % Gaussian kernels
% lddmmoptions.energyweight = [1 1e4]; % weighting between energy of curve and match
lddmmoptions.energyweight = [0 100];
% lddmmoptions.energyweight = lddmmoptions.energyweight/sum(lddmmoptions.energyweight);
lddmmoptions.order = 1;

% image options
% imageoptions.samplesPerPoint = 500;
imageoptions.randomSampling = false;
imageoptions.range = 3*lddmmoptions.scale;

% optimization options
optimoptions.tolFun = 1e-8;
optimoptions.tolX = 1e-10;

% output options
% global globalOptions;
% globalOptions.verbose = true;
clf
visoptions.dim = 2;

% data
lddmmoptions.controlPoints = [40.0; 40.0] % control points
IF = ones(80,80); IF(25:55,25:55) = zeros(31,31);
IM = IF;

options = getDefaultOptions();

% visualize
visoptions.margin = 40;

% pretransform moving image
preScaling = 0.7;
imageoptions.movingTransform = @(x) preScaling*(x-repmat(lddmmoptions.controlPoints,1,size(x,2)))+repmat(lddmmoptions.controlPoints,1,size(x,2));

[methods lddmmoptions imageoptions] = setupImageLDDMM(IM,IF,imageoptions,lddmmoptions,optimoptions);
visualizer = getImageVisualizer(methods.transport,IM,IF,visoptions,imageoptions);
methods.iterationVisualizer = visualizer;
pointVisualizer = getPointVisualizer(methods.transport,lddmmoptions.controlPoints,[],visoptions);

result = runRegister(methods, options);
reshape(result,[],lddmmoptions.L)

visualizer(result);
figure(11), clf
pointVisualizer(result);
for i=1:5
    figure(i), colormap(gray)
end

pause
% pretransform moving image
v = 1/8*pi; 
Rotv = [cos(v) sin(v); -sin(v) cos(v)];
imageoptions.movingTransform = @(x) Rotv*(x-repmat(lddmmoptions.controlPoints,1,size(x,2)))+repmat(lddmmoptions.controlPoints,1,size(x,2));

[methods lddmmoptions imageoptions] = setupImageLDDMM(IM,IF,imageoptions,lddmmoptions,optimoptions);
visualizer = getImageVisualizer(methods.transport,IM,IF,visoptions,imageoptions);
methods.iterationVisualizer = visualizer;
pointVisualizer = getPointVisualizer(methods.transport,lddmmoptions.controlPoints,[],visoptions);

result = runRegister(methods, options);
reshape(result,[],lddmmoptions.L)

visualizer(result);
figure(11), clf
pointVisualizer(result);

