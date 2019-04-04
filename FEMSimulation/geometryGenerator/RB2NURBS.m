function [ nurbs ] = RB2NURBS( rb )
%Convert curves or surface in rB form to nurbs
    nurbs = rb;
    nurbs.form = 'B-NURBS';
    if iscell(rb.knots)
        for i = 1:length(nurbs.knots)
            nurbs.knots{i} = nurbs.knots{i}./max(nurbs.knots{i});
        end
    else
        nurbs.knots = nurbs.knots/max(nurbs.knots(:));
    end
end

