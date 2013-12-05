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

function visualizer = getImageVisualizerL2(transport,IM,IF,visoptions,imageoptions)

margin = imageoptions.margin;
movingTransform = @(x) x;
fixedTransform = @(x) x;
if isfield(imageoptions,'movingTransform')
    movingTransform = imageoptions.movingTransform;
end
if isfield(imageoptions,'fixedTransform')
    fixedTransform = imageoptions.fixedTransform;
end

[IFG,D1IFG,D2IFG] = smoothIG(IF,imageoptions);
[IMG,D1IMG,D2IMG] = smoothIG(IM,imageoptions);

    function stop = imageVisualizer(x,varargin)  
%         x' % output status
        
        % obtain moving image result
        [gridFixed{1} gridFixed{2} gridFixed{3}] = meshgrid(1:size(IF,2),1:size(IF,1),0);

        % transport
        gridMoving = transport(x,gridFixed);
        moved = transport(x,imageoptions.moving);

        % images

        if getOption(visoptions,'showsmooth',true')
            vals = linSampleI(IMG,D1IMG,D2IMG,movingTransform([reshape(gridMoving{1},numel(IF),1) reshape(gridMoving{2},numel(IF),1)]'),imageoptions);
            IMresult = reshape(vals,size(IF));
            vals = linSampleI(IFG,D1IFG,D2IFG,fixedTransform([reshape(gridFixed{1},numel(IF),1) reshape(gridFixed{2},numel(IF),1)]'),imageoptions);
            IFintrp = reshape(vals,size(IF));
            vals = linSampleI(IMG,D1IMG,D2IMG,movingTransform([reshape(gridFixed{1},numel(IF),1) reshape(gridFixed{2},numel(IF),1)]'),imageoptions);
            IMintrp = reshape(vals,size(IF));
        else
            vals = linSampleI(IM,D1IMG,D2IMG,movingTransform([reshape(gridMoving{1},numel(IF),1) reshape(gridMoving{2},numel(IF),1)]'),imageoptions);
            IMresult = reshape(vals,size(IF));
            vals = linSampleI(IF,D1IFG,D2IFG,fixedTransform([reshape(gridFixed{1},numel(IF),1) reshape(gridFixed{2},numel(IF),1)]'),imageoptions);
            IFintrp = reshape(vals,size(IF));
            vals = linSampleI(IM,D1IMG,D2IMG,movingTransform([reshape(gridFixed{1},numel(IF),1) reshape(gridFixed{2},numel(IF),1)]'),imageoptions);
            IMintrp = reshape(vals,size(IF));            
        end
        
        % display
        figure(1), clf, colormap gray, imagesc(IMintrp), colorbar
        hold on, plot(imageoptions.moving(2,:),imageoptions.moving(1,:),'ro');
        hold on, plot(moved(1,:),moved(2,:),'bo');
        figure(2), clf, colormap gray, imagesc(IFintrp), colorbar        
        hold on, plot(imageoptions.moving(1,:),imageoptions.moving(2,:),'ro');
        figure(3), clf, colormap gray, imagesc(IMresult), colorbar
        hold on, plot(imageoptions.moving(1,:),imageoptions.moving(2,:),'bo');
        figure(4), clf, colormap gray, imagesc(IMresult-IFintrp), colorbar
        hold on, plot(imageoptions.moving(1,:),imageoptions.moving(2,:),'ro');
        figure(5), clf, colormap gray, imagesc(abs(IMresult-IMintrp)), colorbar
        
        fprintf('Image matching results:\n')
        fprintf('pre-match 1-norm IM/IF: %f\n',norm(IMintrp(margin:end-margin,margin:end-margin)-IFintrp(margin:end-margin,margin:end-margin),2)/numel(IFG(margin:end-margin,margin:end-margin)));
        fprintf('post-match 1-norm IM/IF: %f\n',norm(IMresult(margin:end-margin,margin:end-margin)-IFintrp(margin:end-margin,margin:end-margin),2)/numel(IFG(margin:end-margin,margin:end-margin)));
        fprintf('relative decrase (%%): %f\n',(1-norm(IMintrp(margin:end-margin,margin:end-margin)-IFintrp(margin:end-margin,margin:end-margin),2)/norm(IMintrp(margin:end-margin,margin:end-margin)-IFintrp(margin:end-margin,margin:end-margin),2))*100);
        fprintf('change 1-norm: %f\n',norm(IMresult-IMintrp,2)/numel(IFG));
        
        stop = false; % dont stop
    end

visualizer = @imageVisualizer;

end