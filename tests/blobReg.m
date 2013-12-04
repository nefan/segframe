% register blob images

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
lddmmoptions.order = 0;

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
lddmmoptions.NpointsX = 4;
lddmmoptions.NpointsY = 4;
IF = zeros(80,80);
center = [40.0; 40.0]; % center of image (for rescaling)
IF(center(1),center(2)) = 255;
smoothscale = 10;
Gfilter = fspecial('gaussian',3*smoothscale,smoothscale);
IF = imfilter(IF,Gfilter,'same'); % smooth
IM = IF;
% scale and range
spacingX = size(IF,1)/lddmmoptions.NpointsX;
spacingY = size(IF,2)/lddmmoptions.NpointsY;
imageoptions.scale = 2*max(spacingX,spacingY)/2;
lddmmoptions.scale = 2.0*imageoptions.scale; % Gaussian kernels

options = getDefaultOptions();

% visualize
visoptions.margin = 20;

% % pretransform moving image
% preScaling = 0.7;
% imageoptions.movingTransform = @(x) preScaling*(x-repmat(center,1,size(x,2)))+repmat(center,1,size(x,2));
% 
% [methods lddmmoptions imageoptions] = setupImageLDDMM(IM,IF,imageoptions,lddmmoptions,optimoptions);
% visualizer = getImageVisualizerL2(methods.transport,IM,IF,visoptions,imageoptions);
% methods.iterationVisualizer = visualizer;
% pointVisualizer = getPointVisualizer(methods.transport,methods.getMoving(),[],visoptions);
% 
% result = runRegister(methods, options);
% reshape(result,[],lddmmoptions.L)
% 
% visualizer(result);
% figure(11), clf
% pointVisualizer(result);
% for i=1:5
%     figure(i), colormap(gray)
% end
% 
% fprintf('done...\n');
% pause

% pretransform moving image
imageoptions.movingTransform = @(x) x;

IM = zeros(80,80);
IM(center(1)+10,center(2)) = 255;
smoothscale = 10;
Gfilter = fspecial('gaussian',3*smoothscale,smoothscale);
IM = imfilter(IM,Gfilter,'same'); % smooth

[methods lddmmoptions imageoptions] = setupImageLDDMM(IM,IF,imageoptions,lddmmoptions,optimoptions);
visualizer = getImageVisualizerL2(methods.transport,IM,IF,visoptions,imageoptions);
methods.iterationVisualizer = visualizer;
pointVisualizer = getPointVisualizer(methods.transport,methods.getMoving(),[],visoptions);

result = runRegister(methods, options);
reshape(result,[],lddmmoptions.L)

visualizer(result);
figure(11), clf
pointVisualizer(result);

fprintf('done...\n');

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

fprintf('done...\n');