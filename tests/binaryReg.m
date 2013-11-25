% register binary images

% random numbers
% s = RandStream.create('mt19937ar','seed',342484);
% s = RandStream.create('mt19937ar','seed',3853873);
% RandStream.setDefaultStream(s);

% image options
% imageoptions.samplesPerPoint = 500;
% imageoptions.randomSampling = false;
imageoptions.measure = 'L2';

% LDDMM options
clear lddmmoptions
% lddmmoptions.energyweight = [1 1e4]; % weighting between energy of curve and match
lddmmoptions.energyweight = [0 10];
% lddmmoptions.energyweight = lddmmoptions.energyweight/sum(lddmmoptions.energyweight);
lddmmoptions.order = 1;

% optimization options
optimoptions.tolFun = 1e-5;
optimoptions.tolX = 1e-5;
optimoptions.numDiff = false;
optimoptions.derivativeCheck = false;

% output options
% global globalOptions;
% globalOptions.verbose = true;
clf
visoptions.dim = 2;

% data
% control points
lddmmoptions.NpointsX = 3;
lddmmoptions.NpointsY = 3;
IF = zeros(80,80); IF(25:55,25:55) = ones(31,31);
center = [40.0; 40.0]; % center of image (for rescaling)
IM = IF;
% scale and range
% lddmmoptions.scale = max(size(IF,1)/lddmmoptions.NpointsX,size(IF,2)/lddmmoptions.NpointsY); % Gaussian kernels
spacingX = size(IF,1)/lddmmoptions.NpointsX;
spacingY = size(IF,2)/lddmmoptions.NpointsY;
imageoptions.scale = max(spacingX,spacingY)/2;
lddmmoptions.scale = imageoptions.scale; % Gaussian kernels

options = getDefaultOptions();

% visualize
visoptions.margin = 40;

% pretransform moving image
preScaling = 0.7;
imageoptions.movingTransform = @(x) preScaling*(x-repmat(center,1,size(x,2)))+repmat(center,1,size(x,2));

[methods lddmmoptions imageoptions] = setupImageLDDMM(IM,IF,imageoptions,lddmmoptions,optimoptions);
visualizer = getImageVisualizerL2(methods.transport,IM,IF,visoptions,imageoptions);
methods.iterationVisualizer = visualizer;
pointVisualizer = getPointVisualizer(methods.transport,methods.getMoving(),[],visoptions);

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
imageoptions.movingTransform = @(x) Rotv*(x-repmat(center,1,size(x,2)))+repmat(center,1,size(x,2));

[methods lddmmoptions imageoptions] = setupImageLDDMM(IM,IF,imageoptions,lddmmoptions,optimoptions);
visualizer = getImageVisualizerL2(methods.transport,IM,IF,visoptions,imageoptions);
methods.iterationVisualizer = visualizer;
pointVisualizer = getPointVisualizer(methods.transport,methods.getMoving(),[],visoptions);

result = runRegister(methods, options);
reshape(result,[],lddmmoptions.L)

visualizer(result);
figure(11), clf
pointVisualizer(result);

