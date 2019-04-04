function [ nurb_surface ] = rb2surface( rb1, rb2, show )
%Rule surface from two rational B splines
%   Convert splines to nurbs first, then rule nurb surface
%   Need Nurbs toolbox
%   http://www.mathworks.com/matlabcentral/fileexchange/26390-nurbs-toolbox-by-d-m-spink
if ~exist('show', 'var') || isempty(show)
    show = false; end

    nurb1 = RB2NURBS(rb1);
    nurb2 = RB2NURBS(rb2);
    nurb_surface = nrbruled(nurb1,nurb2);
    if show == true
        nrbplot(nurb_surface, [100,100]);
        hold on;
    end
    
end

