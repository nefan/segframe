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

function v = gradEnergy(ytt,rhott,L,R,scales,scaleweight,Ks,D1Ks,D2Ks)
%
% Compute energy gradient
% Only Gaussian kernels supported in optimized C function
%

global dim;
CSP = dim*(1+R); % column size particles
CSD = dim*R; % column size directions

% caches
% reverting from containers.Map because matlab2007 doesn't support it
% ytMap = containers.Map('KeyType','double','ValueType','any');
% rhotMap = containers.Map('KeyType','double','ValueType','any');

% counters
lookups = 0;
misses = 0;

% functions below constitute primitive cache since early matlabs doesn't
% have containers.Map
ytind = [];
ytc = cell(1);
rhotind = [];
rhotc = cell(1);
function [i i1] = binsearch(t,c)
    i1 = 0;
    i2 = size(c,1);
    if i2 == 0
        i = -1;
        return;
    end
    vv = c(i2,1);
    if vv == t
        i = i2;
        return;
    else if vv < t
            i = -1;
            i1 = i2;
            return;
        end
    end
    while true
        in = i1+ceil((i2-i1)/2);
        vv = c(in,1);
        if vv == t
            i = in;
            return;
        else if vv < t
                i1 = in;
            else % vv < t
                i2 = in;
            end
        end
        
        if ~(i1 < i2-1)
            break;
        end
    end
    i = -1;
    return; % not found
end
function i = ytcacheddeval(t)
    [i insert] = binsearch(t,ytind);
    lookups = lookups + 1;
    if i ~= -1 % found in cache
        i = ytind(i,2);
        return;
    else % missing
        misses = misses + 1;
        i = size(ytc,2)+1;
        ytind = [ytind(1:insert,:); [t i]; ytind(insert+1:end,:)];
        ytc{i} = deval(ytt,t);
    end
end
function i = rhotcacheddeval(t)
    [i insert] = binsearch(t,rhotind);
    lookups = lookups + 1;
    if i ~= -1 % found in cache
        i = rhotind(i,2);
        return;
    else % missing
        misses = misses + 1;
        i = size(rhotc,2)+1;
        rhotind = [rhotind(1:insert,:); [t i]; rhotind(insert+1:end,:)];
        rhotc{i} = deval(rhott,t);
    end
end

global testC;
function w = Egradtc(tt) % wrapper for C version of Egradt
    w = zeros(size(tt));

    for k = 1:length(tt)
        t = tt(k);
        
        % cache stuff
%         lookups = lookups + 1;
%         if isKey(ytMap,t)
%             yt = ytMap(t);
%         else
%             yt = deval(ytt,t);
%             ytMap(t) = yt;
%             misses = misses + 1;
%         end
         yt = ytc{ytcacheddeval(t)};
%         assert(norm(yt2-yt) < 10e-15);
%         lookups = lookups + 1;
%         if isKey(rhotMap,t)
%             rhot = rhotMap(t);
%         else
%             rhot = deval(rhott,t);
%             rhotMap(t) = rhot;
%             misses = misses + 1;
%         end
         rhot = rhotc{rhotcacheddeval(t)};
%          assert(norm(rhot2-rhot) < 10e-15);

        w(k) = gradE(t,yt,rhot,L,R,dim,scales.^2,scaleweight.^2,GL,GS,GI);
    end
    
    % debug
    if testC
        w2 = Egradt(tt);
        assert(norm(w-w2) < 10e-12);
    end
end

% integrate energy gradient, slow version
function w = Egradt(tt)
    yttt = deval(ytt,tt);

    w = zeros(size(tt));

    for k = 1:length(tt)
        t = tt(k);
        yt = reshape(yttt(:,k),CSD*L,CSP*L);
        rhot = reshape(deval(rhott,t),CSP,L);
    
        for i = 1:L % particle
            xi = rhot(1:dim,i);
            dxi = yt(CSD*(GL-1)+dim*(GS-1)+GI,CSP*(i-1)+(1:dim))';
            
            for ll = 1:L % particle
                xl = rhot(1:dim,ll);
                dxl = yt(CSD*(GL-1)+dim*(GS-1)+GI,CSP*(ll-1)+(1:dim))';
                
                for ss = 1:R                                           
                    al = rhot(dim*(ss-1)+(1+dim:2*dim),ll); 
                    ai = rhot(dim*(ss-1)+(1+dim:2*dim),i);                    
                                                                                    
                    dai = yt(CSD*(GL-1)+dim*(GS-1)+GI,CSP*(i-1)+dim*(ss-1)+(1+dim:2*dim))';
                    dal = yt(CSD*(GL-1)+dim*(GS-1)+GI,CSP*(ll-1)+dim*(ss-1)+(1+dim:2*dim))';                    
    
                    ks = Ks(xi,xl,scales(ss),scaleweight(ss));
                    
                    % de
                    w(k) = w(k) ...
                        + (ks*al)'*dai ...
                        + ai'*(ks*dal ...
                        + 2*al*D1Ks(xi,xl,scales(ss),scaleweight(ss))*(xi-xl)'*(dxi-dxl));
                end
            end
        end
    end
end

% main loop
v = zeros(dim*R*L,1);
for l = 1:L % particles
    for s = 1:R % scales
        GL = l;
        GS = s;
        for GI = 1:dim
            v(dim*R*(l-1)+dim*(s-1)+GI,1) = quad(@Egradtc,0,1); % C version
            %v(dim*R*(l-1)+dim*(s-1)+GI,1) = quad(@Egradt,0,1); % matlab version
        end
    end
end

global verbose
if verbose
    misses
    lookups
end

end