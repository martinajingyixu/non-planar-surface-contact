function [ model ] = igesmodel( R1s, R2s, heights, show ,rshift)
%Build model (circles or ellipses for each height) based on splines for IGES files
%   R1s: vector of ellipse radius (long)
%   R2s: vector of ellipse radius (short)
%   heights: height of the ellipse
if ~exist('show', 'var') || isempty(show)
    show = false; end

if ~exist('rshift', 'var') || isempty(rshift)
    rshift = zeros(size(R1s)); end
% show = true;
    numHeights = length(heights);
    angles = {[0,pi/2],[pi/2,pi],[pi,3*pi/2],[3*pi/2,2*pi]};
    arcs = cell(4,numHeights);
    shifted_arcs = cell(4,numHeights);
    
    if show == true
        figure;
    end
%     figure(1)
% an object consists of multiple ellipse
% create 4 arcs for an ellipse

    for n = 1:numHeights
        for iArcs = 1:4
            arcs{iArcs,n} = fn2fm(rsmak('arc',1,[0;0],angles{iArcs}),'rB');
            shifted_arcs{iArcs,n} = fncmb(fncmb(arcs{iArcs,n},[R1s(n) 0;0 R2s(n);0 0]),[rshift(n);0;heights(n)]);
            arcs{iArcs,n} = fncmb(fncmb(arcs{iArcs,n},[1 0;0 1;0 0]),[0;0;heights(n)]); 
            if show == true
%                 fnplt(arcs{i,n}); 
%                 hold on;
                fnplt(shifted_arcs{iArcs,n},'black');
                hold on;
                axis equal;
            end
    
        end
    end
%     if show == true
%      savefig([filename '.fig']);
%      close all;
%     end

    model = shifted_arcs;
end

