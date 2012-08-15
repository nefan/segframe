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

% test for point LDDMM registration

% LDDMM options
clear lddmmoptions
lddmmoptions.energyweight = [1 8]; % weighting between energy of curve and match
lddmmoptions.energyweight = lddmmoptions.energyweight/sum(lddmmoptions.energyweight);

% scale, Gaussian kernels
lddmmoptions.scales = [1.5];  
% lddmmoptions.scales = [0.2 1.5 3.0];

% sparsity
lddmmoptions.sparsity = false;
lddmmoptions.sparseoptions.alpha = 0;
% 0.007 for point examples, 0.055 for LDDMM 50, 0.20 for LDDKBM 1,50,90 75 points, 0.015 for LDDKBM 10,50,90 256 points
lddmmoptions.sparseoptions.lambdasp = 0.007; 
lddmmoptions.sparseoptions.c = 0.025;

% output options
% global globalOptions;
% globalOptions.verbose = true;
clear visoptions
visoptions.dim = 2;

% moving points
moving = [4.0 2.0;
    4.2 2.2;                    
    4.4 2.3;                            
    4.6 2.4;            
    4.8 2.5;          
    5.0 2.55;          
    5.2 2.5;          
    5.4 2.4;
    5.6 2.3;
    5.8 2.2;        
    6.0 2.0]';    
% fixed - matching against these
fixed = [4.0 5.0;
     4.2 4.8;
     4.4 4.7;        
     4.6 4.6;
     4.8 4.5;          
     5.0 4.45;          
     5.2 4.5;         
     5.4 4.6;
     5.6 4.7;
     5.8 4.8;         
     6.0 5.0]';

options = getDefaultOptions();

[methods lddmmoptions] = setupPointLDDMM(moving,fixed,[],lddmmoptions);

result = runRegister(methods, options);

% visualize
visualizer = getPointLDDMMVisualizer(methods.transport,moving,fixed,visoptions,lddmmoptions);
visualizer(result);
