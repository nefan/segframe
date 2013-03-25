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

function visualizer = getImageVisualizer(transport,IM,IF,visoptions,imageoptions)

movingTransform = @(x) x;
fixedTransform = @(x) x;
if isfield(imageoptions,'movingTransform')
    movingTransform = imageoptions.movingTransform;
end
if isfield(imageoptions,'fixedTransform')
    fixedTransform = imageoptions.fixedTransform;
end

    function stop = imageVisualizer(x,varargin)  
%         x' % output status
        
        % obtain moving image result
        [gridFixed{1} gridFixed{2} gridFixed{3}] = meshgrid(1:size(IF,2),1:size(IF,1),0);

        % transport
        gridMoving = transport(x,gridFixed);

        % spline
        vals = Spline2D(movingTransform([reshape(gridMoving{2},numel(IF),1) reshape(gridMoving{1},numel(IF),1)]')',...
            zeros(numel(IF),1),IM,[0 0],double([1 1]),1);
        IMresult = reshape(vals,size(IF));
        vals = Spline2D(fixedTransform([reshape(gridFixed{2},numel(IF),1) reshape(gridFixed{1},numel(IF),1)]')',...
            zeros(numel(IF),1),IF,[0 0],double([1 1]),1);
        IFspline = reshape(vals,size(IF));
        vals = Spline2D(movingTransform([reshape(gridFixed{2},numel(IF),1) reshape(gridFixed{1},numel(IF),1)]')',...
            zeros(numel(IF),1),IM,[0 0],double([1 1]),1);
        IMspline = reshape(vals,size(IF));
        
        % display
        figure(1), clf, colormap gray, imagesc(IMspline), colorbar
        figure(2), clf, colormap gray, imagesc(IFspline), colorbar        
        figure(3), clf, colormap gray, imagesc(IMresult), colorbar
        figure(4), clf, colormap gray, imagesc(abs(IMresult-IFspline)), colorbar
        figure(5), clf, colormap gray, imagesc(abs(IMresult-IMspline)), colorbar
%         figure(6), clf, imagesc(abs(IMspline-IFspline))        
        
        fprintf('Image matching results:\n')
        fprintf('pre-match 1-norm IM/IF: %f\n',norm(IMspline-IFspline,1)/numel(IF));
        fprintf('post-match 1-norm IM/IF: %f\n',norm(IMresult-IFspline,1)/numel(IF));
        fprintf('relative decrase (%%): %f\n',(1-norm(IMresult-IFspline,1)/norm(IMspline-IFspline,1))*100);
        fprintf('change 1-norm: %f\n',norm(IMresult-IMspline,1)/numel(IF));
        
        stop = false; % dont stop
    end

visualizer = @imageVisualizer;

end