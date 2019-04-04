function [info] = build_model(type, h, r1, r2, res, thick, parameters, options)
% Create artificial object models (curved bottle) as pointcloud/meshes with export
% 2014-04-02
% units in mm; export in m
% res: resolution/sampling_density
% options.tri_select: Index of selected triangle; if not given, user chooses in GUI
% options.rim: Thicker rim (reinforcement) on the top
% options.flat: Flat object along Y-direction: For any Y > r2, Y = r2
% options.dent: [0...1] Carve a valley-like dent oriented along x-axis on one side at given rel. height
% options.closetop: [0/1] Create closed object (no opening on top)
% options.noshow: Do not show anything (GUI)
% options.filename
% options.dir
% options.material_e:  
% options.material_nu:
% options.material_density:
% options.push_features:      Push only on feature points; [ #along-height, #along-cols]
%                             If not given, push set is all outside Verts

global CONSTANTS;

if nargin == 0
    type = 'bottle'; h = 20; r1 = 4; r2 = 2;  
    thick = 0.2; res = 0.5; options = 'noshow';
end
% Check if isfield and evaluate
fset = @(s,f) (isfield(s,f) && any(s.(f)));

info = struct();
info.err = []; force_val=[];
name=[]; X=[]; Vert_outside=[];
Tri_fixed=[]; Tri_free=[]; Tri_outside=[];  Tri_force=[]; 
name = type;
name_opt = '';
if fset(options, 'rim'); name_opt = [name_opt '-rim']; end
if ~isfield(options, 'filename')
    options.filename = sprintf('%s%s-%dx%dx%d-d%d-r%.0f', type, name_opt, h, r1, r2, round(thick*10), res);
end
if ~isfield(options, 'dir')
    options.dir = 'models/'; end
if ~isfield(options, 'dir_iges')
    options.dir_iges = [options.dir '/iges_models/']; end
name_export = [ options.dir '/' options.filename ];
createdir(options.dir_iges);

Tri_force = [];
force_val = [];

% load shapes % Spline 1 calculated with fit; set coefs to 0 for straight lines
X = [];     % 3D points; each row, x y z
rows = [];  % 1xrows(X) vector. Groups equal "height rows" in X for meshing
sample_disk_single = fset(options, 'flat'); % Workaround, see flat-test below 

switch type
    case 'bottle-gen'
        % Generate revolution curve
        % Parameters: see line PARAMETERS below and XXX.svg
        tmp = num2cell(parameters);
        % From bottom to top. H<n>/R<n> are up to 4 convex (R>0) or concave (R<0)
        % rings of height H<n>
        % Neck height and radius
        % Tapering (Verjuengung) to neck has autoheight ancan be 
        %     linear (Rtap=0), convex (Rtap>0), concave (Rtap<0)
        % Shift is a sidewards spline displacement along y-axis
        [riH(1), riR(1),  riH(2), riR(2),  riH(3), riR(3),  riH(4), riR(4), ...
            Hneck, Rneck, Rtap, Rshift] = tmp{:}; % PARAMETERS
        s = cell(0); 
        curbreak = 0;
        % Bottle is made of several rings with splines as rev. curves
        % Iterate over rings:
        for i = 1:length(riH)
            if riH(i) == 0;
                continue; end
            s{end+1} = spline(curbreak + [0 riH(i)/2 riH(i)], [1 1+riR(i) 1]);
            curbreak = s{end}.breaks(end);
        end
        % Tapering (Verjuengung)
        Htmp = [curbreak, nan, 1-Hneck]; Htmp(2) = mean(Htmp([1,3]));
        if Htmp(1) < Htmp(3) && Hneck > 0
            s{end+1} = spline(Htmp, [1 (1+Rneck)/2+Rtap Rneck]); end
        % Neck
        if Hneck > 0
            s{end+1} = spline([1-Hneck 1], [Rneck Rneck]); end
        % Sideward shift
        sshift = spline([0 0.5 1], [0 Rshift 0]);
        
        % Put the splines together
        sall = struct('form', 'pp', 'breaks', s{1}.breaks(1), ...
            'coefs', zeros(0,3), 'order', 3, 'dim', 1, 'pieces', 0);
        for si = 1:length(s)
            sc = s{si};
            if sall.breaks(end) ~=  sc.breaks(1); warning('bottle-gen: Bad spline'); end
            sall.breaks = [ sall.breaks(1:end-1) sc.breaks ];
            sall.coefs(end+1, : ) = 0;
            sall.coefs(end, end-size(sc.coefs,2)+1:end)=  sc.coefs;
            sall.pieces = sall.pieces + sc.pieces;
        end
        if fset(options, 'flat');
            r3 = r1; else r3 = r2; end
        
        if fset(options, 'rim')
            heights_iges = 0:res:(h+res);
        else
            heights_iges = 0:res:h;
        end
        r1_iges_ = r1 * ppval(sall, heights_iges/h);
        r2_iges = r3 * ppval(sall, heights_iges/h);
        rshift = r1 * ppval(sshift, heights_iges/h);
        
        
        rbsplines = igesmodel([r1_iges_],[r2_iges],[heights_iges],~options.noshow,rshift);
        squeezeLocationArray.xaxis = [[r1_iges_];zeros(1,length([r1_iges_]));[heights_iges]];
        squeezeLocationArray.yaxis = [zeros(1,length([r2_iges]));[r2_iges];[heights_iges]];
%         hold on
        model_ruledSurfaces = rulesurf(rbsplines, ~options.noshow);
%         bottemSurface = ruleEllipse( r1,r2,0,~options.noshow );
        %bottemSurface = ruleEllipse( r1_(1),r2_(1),0,~options.noshow );
        %model_ruledSurfaces{end+1} = bottemSurface;
        model_splines = rbsplines(:);
        if CONSTANTS.SAVEIMAGES
        print(gcf, ['/home/jingyi/promotion/programming/softContactAnalysis/dataset/figures/fig' num2str(CONSTANTS.FILENAME) '.pdf'], '-dpdf','-bestfit');
        CONSTANTS.FILENAME = CONSTANTS.FILENAME+1;
        end
        X = [0 0 0]; rows = [1]; % Center of bottom plate
        heights = res:res:h;
        r1_ = r1 * ppval(sall, heights/h);
        r2_ = r3 * ppval(sall, heights/h);
        yshift = r1 * ppval(sshift, heights/h);

        [X_, rows_] = sample_disk(r1, r3, res, sample_disk_single); % Bottom plate
        X = [ X ; X_ ]; rows = [rows rows_+1];
        for i = 1:length(heights)
            X_ = sample_ellipse(r1_(i), r2_(i), res);
            X_(:,3) = heights(i);
            X_(:,1) = X_(:,1) + yshift(i);
            rows = [ rows , (rows(end)+1) * ones(1, size(X_,1)) ];
            X = [ X ; X_ ];
        end
        

    case 'bowl'
        % r1/r2:   Axes of elliptic base
        % parameters(1): Radius of bottom plate, relative to r1
        % parameters(2): Radius at half height, realative. =0 linear, >0 concave, <0 convex
        botR = parameters(1);
        halfR = parameters(2);
        X = [ 0 0 0 ];
        rows = [ 1 ];
        heights = 0:res:h;
        sall = spline([0 0.5 1], [botR  (botR+1)/2+halfR  1]);
        heights_iges = heights;
        if fset(options, 'rim')
            heights_iges = 0:res:(h+res);
        end
        % Revolution curve
        % old: r = r1 * (botR + (1-botR) * (heights/h)) .^ (0.5);
        r1_iges = r1 * ppval(sall, heights_iges/h);
        r2_iges = r2 * ppval(sall, heights_iges/h);

        
        % Evaluate splines, etc.
        X = [0 0 0]; rows = [1]; % Center of bottom plate
        rshift = zeros(size(r1_iges));
        sshift = 0;
        rbsplines = igesmodel([r1_iges],[r2_iges],[heights_iges],~options.noshow);
%         hold on
        model_ruledSurfaces = rulesurf(rbsplines, ~options.noshow);
        
        squeezeLocationArray.xaxis = [[r1_iges];zeros(1,length([r2_iges]));[heights_iges]];
        squeezeLocationArray.yaxis = [zeros(1,length([r1_iges]));[r2_iges];[heights_iges]];
        if CONSTANTS.SAVEIMAGES
        print(gcf, ['/home/jingyi/promotion/programming/softContactAnalysis/dataset/figures/fig' num2str(CONSTANTS.FILENAME) '.pdf'], '-dpdf','-bestfit');
        CONSTANTS.FILENAME = CONSTANTS.FILENAME+1;
        end

        %bottemSurface = ruleEllipse( r(1),r(1),0,~options.noshow );
        %model_ruledSurfaces{end+1} = bottemSurface;
        model_splines = rbsplines(:);
        
        r = r1 * ppval(sall, heights/h);
        for i = 1:length(heights)
            X_ = sample_ellipse(r(i), r2/r1*r(i), res);
            X_(:,3) = heights(i);
            rows = [ rows , (rows(end)+1) * ones(1, size(X_,1)) ];
            X = [ X ; X_ ];
        end
        

end
% ------ Other geometry-related options ------
% Flat object
if fset(options, 'flat')
    % Special case disk on floor: It has multiple points at same height which
    % would be projected onto the same line y = +-r2.
    % Workaround: sample_disk (called above) creates only 1 ellipse.
    % Todo: check for all models!
    X(X(:,2) > r2, 2) = r2;
    X(X(:,2) < -r2, 2) = -r2;
end

% Dent in object (like a vallay along x-axis)
% options.dent gives relative height, rest is hard-coded here
if fset(options, 'dent')
    % Depth of dent
    d_dent = r2 * 0.2;
    % Check at which height to create dent
    h_dent = options.dent * h;
    [~, hi_dent] = min(abs(h_dent - heights));
    % Find vertices at this height and max y-value:
    mask_h = (X(:,3) == heights(hi_dent));
    ymax = max(X(mask_h,2));
    % Valley height of dent:
    y_dent = ymax - d_dent;
    mask_hy = (mask_h & (X(:,2) > y_dent)) > 0;
    if sum(mask_hy) == 0;
        warning('Found no vertices for dent!'); end
    % Carve the dent
    X(mask_hy, 2) = y_dent;
end

% Add rim on top; this is maybe obsolete with stiff rim-material in Vega
mr = max(rows);
toprim = rows == mr;
if fset(options, 'rim');
    X_ = X(toprim,:);
    % Thickness of rim; limit by smallest radius
    thick2 = min([ thick*10, min(max(abs(X_(:,1:2))))/2 ]);
    X_(:,3) = X_(:,3) + res;
    X = [ X ; X_ ];
    rows(1, end+1:end+size(X_,1)) = max(rows) + 1;
    toprim = rows >= mr;
else
    thick2 = thick;
end

% Close the opening on the top
if fset(options, 'closetop');
    X = [ X ; 0 0 X(end,3) ];
    rows = [ rows , (rows(end)+1)*[1] ];
end

% ----- Create triangle mesh -----
[T, TriN] = mesh_structured(X, rows); % T: Reference (index) to 3 vertices for 1 triangle
%save([name_export, '_manuell_sampling.mat'],'X','rows');

% Find the floor trianges (z=0) for Dirichlet condition
% JY: If the z-coords. of the vertix of a triangle is 0, then they are on
% the floor and are fixed.
i_floor = 0;
for i = 1:size(T,1)
    if any(X(T(i,:), 3) ~= 0)
        break
    end
    i_floor = i;
end
assert(i_floor > 0);
Tri_fixed = T(1:i_floor, :); % the fixed triangles on the ground


% ----- Inner/outer surface -----
Xo = X;  % o: outer surface
No = find_normals(Xo, T, TriN);
% if ~fset(options, 'noshow')
%     figure(1); clf; hold on;
%     l = res/2;
%     tr = @(x) reshape(x', 1,[]);
%     no_plot = [ tr(Xo);  tr(Xo + l*No) ];
% %     plot3(no_plot(:,1:3:end), no_plot(:,2:3:end), no_plot(:,3:3:end), ...
% %         'b', 'LineWidth', 3);
%     axis equal
% end
%Xi = inner_surface(X(~toprim,:), ctr, thick);
%Xit = inner_surface(X(toprim,:), ctr, thick2);
if thick > 0
    Xi  = X - thick * No;  % i: inner surface
    Xi(toprim,:) = X(toprim,:) - thick2 * No(toprim,  :);
    % Xit = X(toprim,:)  - thick2 * No(toprim,  :);
    X = [ Xo ; Xi ];
    assert(size(X,1) == 2*size(Xo,1));
    % Mesh inner surface just like outer surface, but invert triangle
    % orientation
    Ti = fliplr(T) + size(Xo, 1);
    if ~fset(options, 'closetop');
        % Mesh the gap btwn innner and outer surface seen from top
        toprow = rows == max(rows);
        rows_gap = [ toprow , toprow*2 ];
        Tgap = mesh_structured(X, rows_gap);
    else Tgap = [];
    end
    Tall = [ T ; Ti ; Tgap];
    Nall = [ No ; -No ];
else
    Nall = No;
    Tall = T;
end

% Check for NANs
nnan = sum(isnan(X(:)));
if nnan > 0;
    info.err= [info.err, sprintf('vertex-nans:%d,', nnan)]; end
% 
% % ----- Plot; Select 1 triangle via GUI -----
% if ~fset(options, 'noshow')
%     figure(1);
%     xlabel X, ylabel Y, zlabel Z
%     colormap(hsv(32));
%     Cidx = mod(1:size(T,1), 32) + 1;
%     %C = col(idx, :);
%     h_tri = zeros(1, size(T, 1));
%     for i = 1:size(T, 1)
%         h_tri(i) = trisurf(T(i,:), X(:,1), X(:,2), X(:,3), ...
%             'HitTest', 'on', 'FaceAlpha', 0.9);
%     end
%     %trisurf(T, X(:,1), X(:,2), X(:,3), Cidx, 'HitTest', 'on');
%     
% %     plot3(Xo(:,1), Xo(:,2), Xo(:,3), 'r.');
% % 
% %     if exist('Xi', 'var')
% %         plot3(Xi(:,1), Xi(:,2), Xi(:,3), 'g.'); end
%     view(45,60);
%     axis equal
% end

if ~isfield(options, 'tri_select')
    if ~fset(options, 'noshow')
        axis equal
        title('Select one triangle by clicking and press enter');
        pause
        Tsel = find(gco == h_tri);
        Tsel, 
        trisurf(T(Tsel,:), X(:,1), X(:,2), X(:,3), Cidx(Tsel), ...
            'EdgeColor', [1 0 0], 'FaceColor', [0 0 0]);
        title('');
    else
        Tsel = 0;
    end
else
    Tsel = options.tri_select;
end

% Selected triangle
isel = i_floor+1:size(Tall,1);
if Tsel > 0
    Tri_force = T(Tsel, :);
    force_val = -No(T(Tsel,1), :);
    assert(any(isel == Tsel));
    isel(isel == Tsel) = [];
end

% ----- Different sets of triangles/vertices, feature locations -----
% Free triangles: Neumann-surface, no force
Tri_free = Tall(isel, :);
Tri_outside = T(i_floor+1:end, :);
fixed_vertics = unique(Tri_fixed);
tmp_free_vertics = unique(Tri_outside);
Vert_outside = setdiff(tmp_free_vertics,fixed_vertics);
if fset(options, 'push_features')
    feat_idx = feat_locations(X(Vert_outside, :), options.push_features(1), options.push_features(2));
    Vert_push = unique(Vert_outside(feat_idx));
else
    Vert_push = Vert_outside;
end

if ~fset(options, 'skiptet')

end %if ~fset(options, 'skiptet')

% ----- Export -----
scale = 0.001; % Should probably be 0.001 for meters
fprintf('Export to %s [.igs,.obj,.mat]\n', name_export);
models2iges = [model_ruledSurfaces(:);model_splines(:)];
igesout(models2iges, options.dir_iges, options.filename);

% 3D surface mesh obj (also used by Vega), [stl]
fh = fopen([name_export '.obj'], 'w');
fprintf(fh, 'v %f %f %f\n', scale*X');
fprintf(fh, 'f %d %d %d\n', [ Tri_fixed ; Tri_free ]');
fclose(fh);

if thick > 0
    save([name_export '.mat'], 'X', 'Tri_fixed', 'Ti', 'Tri_free', 'options', ...
    'Tri_outside', 'Tri_force', 'Vert_outside', 'Vert_push', 'No', 'TriN','squeezeLocationArray','rshift','sall','sshift');
else
    save([name_export '.mat'], 'X', 'Tri_fixed', 'Tri_free', 'options', ...
    'Tri_outside', 'Tri_force', 'Vert_outside', 'Vert_push', 'No', 'TriN','squeezeLocationArray','rshift','sall','sshift');
end
if ~isempty(info.err)
    warning(['Errors occured: ' info.err]); end
return


function [T,n] = mesh_structured(X, rows)
% Meshing with triangles, of structured points, row-by-row.
% A standard approach - like ball pivoting, etc - might be a better solution
% Each point i in current row spans triangle:
%    [ Current point, Next point in current row, closest point in previous row ]
lasti = [];
T = [];  % point indexes for triangles
urows = unique(rows); % Unique height rows
urows(urows==0) = [];
for ri = urows
    %if ri == urows(end); keyboard; end
    cur = rows == ri;
    curi = find(cur); % All points at current height level
    curi(end+1) = curi(1);
    last_besti = []; first_besti = [];
    if ~isempty(lasti)
        for ic = 1:length(curi)-1 % Over all points in current height row
            ix = curi(ic);
            ix_right = curi(ic+1);
            xc = mean(X([ix ix_right],:),1); % Point btwn vertices [ix ix_right]
            % find besti: closest point btwn X(lasti,:) and xc
            d = sum((repmat(xc, length(lasti), 1) - X(lasti, :)) .^ 2, 2);
            [~,di] = min(d);
            assert(di > 0);
            besti = lasti(di);
            % Add triangle: ix - ix_right - besti
            if ix ~= ix_right % normal case
                T(end+1,:) = [ ix ix_right besti ];
            else % Special case: current height row has only one point
                lasti_ = [lasti lasti(1)];
                for il = 1:numel(lasti_)-1
                    T(end+1,:) = [ ix lasti_(il) lasti_(il+1) ];
                end
            end
            % Add more triangles with points ix and anything btwn [last_besti..besti]
            if ~isempty(last_besti)
                bii = find(last_besti == lasti);
                while lasti(bii) ~= besti
                    bii_ = bii + 1;
                    if bii_ > length(lasti); bii_ = 1; end
                    T(end+1,:) = [ ix lasti(bii) lasti(bii_) ];
                    bii = bii_;
                    % if lasti(bii) == besti; break; end
                end
            end
            % If last loop: Another triangle to fill last hole
            % Tringle is [ ix_right(cur row)  besti (last row)  first_besti (last row) ]
            % if no points btwn besti and first_besti
            if ix == curi(end-1) && ~isempty(first_besti) && first_besti ~= besti
                bii = find(besti == lasti);
                % Loop over all points [ besti ... first_besti ] in last row
                while 1;  % lasti(bii) ~= first_besti
                    bii_ = bii + 1;
                    if bii_ > length(lasti); bii_ = 1; end
                    T(end+1,:) = [ ix_right lasti(bii) lasti(bii_) ];
                    bii = bii_;
                    if lasti(bii) == first_besti; break; end
                end
            end
            %
            if ix == curi(1)
                first_besti = besti;
            end
            last_besti = besti;
        end
    end
    lasti = curi(1:end-1); % last height row
end
% Reorient triangles for consistent orientation
% Calculate triangle normals
n = cross(X(T(:,2), :) - X(T(:,1), :), X(T(:,3), :) - X(T(:,1),:));
nl = sqrt(sum(n.^2, 2));
n = n ./ [nl, nl, nl];
% Check if normals point inward in xy-plane
X_ = X(T(:,1), :);
Xn_ = X_ + 0.001 * n;
X_(:,3) = 0; Xn_(:,3) = 0;
n_inside = sum(Xn_.^2, 2) < sum(X_.^2, 2);
% For triangles oriented (mostly) along x-y plane:
X_ = X(T(:,1), :);
n_xy = abs(n(:,3)) > cos(30/180*pi);
n_inside(n_xy) = xor(n(n_xy, 3) > 0, (X_(n_xy, 3) > 0.01));
% Flip direction of triangle
T(n_inside, :) = fliplr(T(n_inside, :));
n(n_inside, :) = -n(n_inside, :);
return

function [X, rows] = sample_disk(a, b, res, single)
if ~exist('single', 'var')
    single = 0; end
X = [];
rows = [];
if ~single
    d_ = min(a,b)-res:-res:0; end
if single || isempty(d_)
    d_ = 0; end
for di = 1:length(d_)
    d = d_(di);
    X_ = sample_ellipse(a-d, b-d, res);
    if size(X_,1) <= 1; continue; end
    X = [ X ; X_ ];
    rows = [ rows, repmat(di, 1, size(X_,1)) ];
end
%plot(X(:,1), X(:,2), 'r.');
%axis equal;
return


function [X, rows] = sample_disk_rect(r1, r2, res)
%sample the bottom and up part

    X = []; rows = [0];

    r1s = -r1:res:r1;
    r2s = -r2:res:r2;

    for r1i = 1:length(r1s)
        X_ = [];
        for r2i = 1:length(r2s)                    
            r1_ = r1s(r1i);
            r2_ = r2s(r2i);
            X_ = [X_;r1_ r2_ 0]; % group all X_ with same height together 
        end
        rows = [ rows , (rows(end)+1) * ones(1, size(X_,1)) ];
        X = [ X ; X_ ];
    end    
    rows = rows(2:end);

return

function X = sample_ellipse(a, b, res)
% Find local density of points on ellipse with fixed delta t_:
t_ = linspace(0, 2*pi, 101)'; %generate 101 points between 0 and 2*pi
X_ = diff([ a*cos(t_), b*sin(t_) ]); %diff of the matrix is the matrix of row differences.
dl = sqrt(sum(X_.^2, 2)) * 100;
smpl = dl/res; % Unit sampling density. 1 = 1 sample/2pi
smpl(smpl < 8) = 8;
% smpl(smpl < 25) = 25; %JY: more points
t = [0];
% Find t-values to sample ellipse:
while 1
    i = min(100, 1 + floor(100 * t(end)/(2*pi)));
    t_add = t(end) + 2*pi/smpl(i);
    if t_add >= 2*pi; break
    else; t(end+1) = t_add; end
end
% Sample ellipse:
X = [ a*cos(t'), b*sin(t'), zeros(size(t')) ];
% Start/end too close?
if norm(X(1,:) - X(end,:)) < res/2
    X(end,:) = [];
end
%plot(X(:,1), X(:,2), 'r.'); axis equal;

return



function model_ruledSurfaces = rulesurf(rbsplines, show)
    %rule surface between rb lines. 
    % rbsplines is 4*n cell, each column is an ellipse, consists of 4 arcs.
    [rows,cols] = size(rbsplines);
    N = (rows-1)*cols;
    model_ruledSurfaces = cell(N ,1);
    n = 1;
    if show == 1
%         figure
        axis equal
        hold on
        view(25,30)
    end
    for c = 1:cols-1
        for r = 1:rows
            rb1 = rbsplines(r,c);
            rb2 = rbsplines(r,c+1);
            model_ruledSurfaces{n} =  rb2surface( rb1{1}, rb2{1}, show);
            n = n + 1;
        end
    end
     if show == 1
    hold on
    ax = gca;
%     labelFontSz = 14;
%     ax.FontSize = labelFontSz;
%     xlabel('x','FontSize',labelFontSz)
%     ylabel('y','FontSize',labelFontSz)
%     zlabel('z','FontSize',labelFontSz)
%     grid on
    ax.XRuler.Axle.LineStyle = 'none';  
    axis off

     end
    
    
return


function N = find_normals(X, T, Tn)
N = zeros(size(X));
for i = 1:size(X, 1)
    iT = any(T == i, 2);  % triangles with point i
    n  = sum(Tn(iT, :), 1);
    N(i,:) = n / norm(n);
end
return



