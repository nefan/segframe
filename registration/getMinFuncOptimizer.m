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

function optimizer = getMinFuncOptimizer(varargin)
 % Limitied memory BGFS optimizer

if size(varargin,2) > 0
    userOptions = varargin{1};
else
    userOptions = [];
end

options = [];
options.Method = 'lbfgs';  
if getOption(userOptions,'derivativeCheck')
    options.DerivativeCheck = 'on';
else
    options.DerivativeCheck = 'off';   
end

% tolerance
if isfield(userOptions,'tolFun')
%     options = optimset(options,'TolFun',userOptions.tolFun);
    options.progTol = userOptions.tolFun;
end
if isfield(userOptions,'tolX')
%     options = optimset(options,'TolX',userOptions.tolX);
    options.optTol = userOptions.tolX;
end

if isfield(userOptions,'maxIter')
    options.MaxIter = userOptions.maxIter;
else
    options.MaxIter = 1000;
end
if isfield(userOptions,'maxFunEvals')
    options.MaxFunEvals = userOptions.maxFunEvals;
else
    options.MaxFunEvals = 1000;
end


% output options
if getOption(userOptions,'verbose')
    options.Display = 'full';
%     options = optimset(options,'Diagnostics','on');    
else
    options.Display = 'final';
end

    function res = minFuncOptimizer(initialData, F, iterationVisualizer)        
        options.OutputFcn = iterationVisualizer;
        res = minFunc(F,initialData,options);
    end

optimizer = @minFuncOptimizer;

end