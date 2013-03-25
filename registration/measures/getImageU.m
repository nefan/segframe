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

function gradU = getImageU(IM,IF,moving,imageoptions,dim,L)

scale = imageoptions.scale;
range = imageoptions.range;
order = imageoptions.order;
if order == 0
    dimX1 = dim;
else
    dimX1 = dim+dim^2;
end
movingTransform = @(x) x;
fixedTransform = @(x) x;
if isfield(imageoptions,'movingTransform')
    movingTransform = imageoptions.movingTransform;
end
if isfield(imageoptions,'fixedTransform')
    fixedTransform = imageoptions.fixedTransform;
end

% generate sample points
samples = []; % contains samples (columns) for all L points
randomSampling = getOption(imageoptions,'randomSampling',false);
for l = 1:L
    p = moving(:,l);
    
    if randomSampling
        samplesl = zeros(1+dim+1,imageoptions.samplesPerPoint);

        % pseudorandom numbers
        assert(dim == 2);
        samplesl(2:(1+dim),:) = mvnrnd([0; 0],scale^2*eye(dim),imageoptions.samplesPerPoint)';
        samplesl(1+dim+1,:) = 1/imageoptions.samplesPerPoint*ones(1,imageoptions.samplesPerPoint); % weight        
%         figure(2)
%         plot(samplesl(2,:),samplesl(3,:),'+')
%         hold on, plot(p(1),p(2),'r*')
%         axis equal
    else
        samplesl = [];
        if dim == 2
            delta = 2*range/(2*range)*1;
            for x = -range:delta:range
                for y = -range:delta:range
                    d = sqrt(x^2+y^2);
                    if d < range
                        samplesl = [samplesl [0; x; y; exp(-d^2/scale^2)]];
                    end
                end
            end
        else
            assert(dim == 3);
            delta = 2*range/(2*range)*2;
            samplesl = zeros(1+dim+1,length(-range:delta:range)^3);
            c = 0;
            for x = -range:delta:range
                for y = -range:delta:range
                    for z = -range:delta:range
                        d = sqrt(x^2+y^2+z^2);
                        if d < range
                            c = c + 1;
                            samplesl(:,c) = [0; x; y; z; exp(-d^2/scale^2)];
                        end
                    end                    
                end
            end
            samplesl = samplesl(:,1:c);
        end
        samplesl(1+dim+1,:) = samplesl(1+dim+1,:)/sum(samplesl(1+dim+1,:));
%         Nsamples = dim*20+1;
%         disp = scale/5;
% 
%         samplesl = [0; zeros(dim,1)];
% %         for i = 1:(Nsamples-1)/(2*dim)
%         for i = 1:(Nsamples-1)/2
% %             for j = 1:dim
%              for j = 2
%                 ej = zeros(dim,1); ej(j) = 1;
%                 samplesl = [samplesl [0; i*disp*ej]];
%                 samplesl = [samplesl [0; -i*disp*ej]];
%             end
%         end
%        [xx idx] = sort(samplesl(3,2:end)); % sort for plotting
%        samplesl(2:end,2:end) = samplesl(2:end,1+idx);
    end

    % compute values in fixed image
    if dim == 2
        vals = Spline2D(fixedTransform([p(2)+samplesl(3,:)' p(1)+samplesl(2,:)']')',zeros(size(samplesl,2),1),IF,[0 0],double([1 1]),1);
    else
        assert(dim == 3);
        vals = SplineInterpolation(fixedTransform([p(2)+samplesl(3,:); p(1)+samplesl(2,:); p(3)+samplesl(4,:)])',IF,[0 0 0],[1 1 1]);
    end
    samplesl(1,:) = vals;

    % store
    samplesl = reshape(samplesl,[],1);
    samples(:,l) = samplesl;
end


    function [y v] = lgradU(x)
        %
        % gradient, image matching
        %                
   
        x = reshape(x,dimX1,L);
        y = 0;
        v = zeros(dimX1,L);
        for l = 1:L
            p = x(1:dim,l);
            if order == 1
                A = reshape(x(dim+(1:dim^2),l),dim,dim);
            end
            samplesl = reshape(samples(:,l),1+dim+1,[]);
            Nsamples = size(samplesl,2);
            if order == 0
                if dim == 2
                    evalPoints = [p(1)+samplesl(2,:); p(2)+samplesl(3,:)];
                else
                    assert(dim == 3);
                    evalPoints = [p(1)+samplesl(2,:); p(2)+samplesl(3,:);  p(3)+samplesl(4,:)];
                end
            else
                evalPoints = zeros(dim,Nsamples);
                for n = 1:Nsamples
                    evalPoints(:,n) = p+A*samplesl(2:dim+1,n);
                end
            end

            % measure similarity
            if dim == 2
                [valS d2 d1] = PNorm2D(movingTransform([evalPoints(2,:)' evalPoints(1,:)']')',...
                    samplesl(1,:)',IM,[0 0],double([1 1]),1,samplesl(1+dim+1,:));
                ds = [d1'; d2'];
            else
                assert(dim == 3);
                [valS d2 d1 d3] = PNorm(movingTransform([evalPoints(2,:); evalPoints(1,:); evalPoints(3,:)])',...
                    samplesl(1,:)',IM,[0 0 0],[1 1 1],Nsamples*samplesl(1+dim+1,:),1);                
                ds = [d1'; d2'; d3'];
            end

            % value        
            y = y + valS; % debug +sum(reshape(A.^2,[],1));

            % gradient
            g = sum(ds,2);
            v(1:dim,l) = g;
            if order == 1 % gradient Dphi
%                 figure(2),quiver(samplesl(2,:),samplesl(3,:),ds(1,:),ds(2,:)),axis equal
                v(dim+dim*0+(1:dim),l) = ds*samplesl(2,:)';
                v(dim+dim*1+(1:dim),l) = ds*samplesl(3,:)';
                if dim == 3
                    v(dim+dim*2+(1:dim),l) = ds*samplesl(4,:)';
                end
%                 v(dim+(1:dim^2),l) = reshape(2*A,dim^2,1);

%                 % testing
%                 v(dim+dim*0+1,l) = 0;
%                 v(dim+dim*0+2,l) = 0;
%                 v(dim+dim*1+1,l) = 0;
%                 v(dim+dim*1+2,l) = 0;
            end
        end
        
        y = y/L;
        v = reshape(v,dimX1*L,1)/L;
    end

gradU = @lgradU;

end