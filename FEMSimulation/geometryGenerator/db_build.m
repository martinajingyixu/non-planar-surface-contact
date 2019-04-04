function [] = db_build(dbDir, objectParameters, range)

load([ dbDir 'db.mat' ]); %loads Params, defaults
if ~exist('show', 'var'); show = false; end
if ~exist('range', 'var') || isempty(range)
    range = 1:height(objectParameters); end


name_para = arrayfun(@(i) {['p_' num2str(i)]}, 1:12); % Names of the general parameters
build_infos = cell(height(objectParameters), 1);

for it = range
    cp = table2struct(objectParameters(it, :));
    parameters = table2array(objectParameters(it, name_para));
    fprintf('Model %5d/%5d: %s\n', it, height(objectParameters), cp.filename);
    
    cp.noshow = ~show;
%     cp.noshow = 0;

    cp.dir = dbDir;
    cp.tri_select = 0;
    for i = fieldnames(defaults)'
        if ~isfield(cp, i{1})
            cp.(i{1}) = defaults.(i{1});  end; end
    if cp.enclosed > 0
        cp.closetop = 1; end

    info = build_model(cp.type, cp.h, cp.r1, cp.r2, cp.res, cp.thick, parameters, cp);
    build_infos{it} = info;
    if show; pause; end

end
save([ dbDir '/db.mat' ], 'build_infos', 'objectParameters', '-append');