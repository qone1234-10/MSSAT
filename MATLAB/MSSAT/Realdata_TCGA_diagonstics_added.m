%% =========================================================
%% TCGA MSSAT partial joint projection diagnostics
%% Targets
%%   subset 14: blocks 1 and 4
%%   subset 24: blocks 2 and 4
%% ---------- user paths ----------
if ~exist('save_dir', 'var') || isempty(save_dir)
    save_dir = '~\Matlab\MSSAT\saved_results';
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

if ~exist('save_file', 'var') || isempty(save_file)
    save_file = fullfile(save_dir, 'tcga_4block_analysis.mat');
end

%% ---------- load existing save_struct if needed ----------
if exist(save_file, 'file') == 2
    S = load(save_file, 'save_struct');
    save_struct = S.save_struct;
else
    error('save_file does not exist: %s', save_file);
end

if ~isfield(save_struct, 'realdata')
    error('save_struct.realdata is missing. Run the TCGA loading section first.');
end

realdata = save_struct.realdata;

if isfield(save_struct, 'fields')
    fields = save_struct.fields;
else
    fields = arrayfun(@(k) sprintf('block%d', k), 1:numel(realdata), 'UniformOutput', false);
end

if ~exist('time', 'var')
    if isfield(save_struct, 'time')
        time = save_struct.time;
    else
        time = struct();
    end
end

if ~exist('jrank', 'var')
    if isfield(save_struct, 'jrank')
        jrank = save_struct.jrank;
    else
        jrank = struct();
    end
end

if ~exist('result_real', 'var') && isfield(save_struct, 'result_real')
    result_real = save_struct.result_real;
end

%% ---------- if already computed, just display and skip partial MSSAT rerun ----------
skip_partial_mssat_rerun = false;

if isfield(save_struct, 'mssat') && isstruct(save_struct.mssat) && ...
        isfield(save_struct.mssat, 'partial_projection_table')
    disp('Existing MSSAT partial projection diagnostics found. No mssat rerun.');
    disp(save_struct.mssat.partial_projection_table);

    if isfield(save_struct.mssat, 'partial_projection_angle_table')
        disp('===== Principal angle long table =====');
        disp(save_struct.mssat.partial_projection_angle_table);
    end

    skip_partial_mssat_rerun = true;
end

if ~skip_partial_mssat_rerun

%% =========================================================
%% Sequential partial MSSAT + diagnostics in the same loop
%% =========================================================
disp('================ MSSAT partial + projection diagnostics ================');

X_list = realdata;
K = numel(X_list);
fprintf('Running sequential partially joint mssat on %d data blocks (BRCA)...\n', K);

%% row-centered residual initialization
R_list = cell(1, K);
for k = 1:K
    Xk = double(X_list{k});
    R_list{k} = bsxfun(@minus, Xk, mean(Xk, 2));
end

n  = size(R_list{1}, 2);
In = eye(n);

subset_list = { [1 2 3], [1 2 4], [1 3 4], [2 3 4], ...
                [1 2],   [1 3],   [1 4],   [2 3],   [2 4],   [3 4] };
subset_labels = {'123','124','134','234','12','13','14','23','24','34'};

partial_rank_names = {'r123','r124','r134','r234','r12','r13','r14','r23','r24','r34'};
partial_ranks_vec = zeros(1, numel(subset_list));

%% Save outputs so diagnostics can be reused later without rerun.
partial_outputs = cell(1, numel(subset_list));
partial_P_joint = cell(1, numel(subset_list));
partial_V_joint = cell(1, numel(subset_list));

%% Target diagnostics.
target_labels = {'14', '24'};
target_mask = ismember(subset_labels, target_labels);
partial_projection_diag = struct();

%% Table containers.
diag_subset = {};
diag_block_id = [];
diag_block_name = {};
diag_joint_rank = [];
diag_initial_rank = [];
diag_subspace_projection_error = [];
diag_matrix_projection_error = [];
diag_mean_angle_deg = [];
diag_max_angle_deg = [];
diag_principal_angles_deg = {};

angle_subset = {};
angle_block_id = [];
angle_block_name = {};
angle_index = [];
angle_value_deg = [];

fprintf('\n===== Sequential partially joint rank estimation and target diagnostics =====\n');
tic;

for s = 1:numel(subset_list)
    idx = subset_list{s};
    lab = subset_labels{s};

    X_sub = R_list(idx);

    out_sub = mssat( ...
        X_sub, ...
        'center', true, ...
        'alpha', normcdf(-6), ...
        'rank_method', 'bema', ...
        'sigma_method', 'bema', ...
        'test_mode', 'alpha', ...
        'mode', 'test', ...
        'alt_angle_deg', 10, ...
        'verbose', false);

    r_hat_s = double(out_sub.joint_rank);
    partial_ranks_vec(s) = r_hat_s;

    fprintf('Subset %s: estimated joint rank = %d\n', lab, r_hat_s);

    if isfield(out_sub, 'P_joint') && ~isempty(out_sub.P_joint)
        P_sub = out_sub.P_joint;
    else
        P_sub = zeros(n, n);
    end

    if isfield(out_sub, 'V_joint') && ~isempty(out_sub.V_joint)
        V_sub = out_sub.V_joint;
    else
        V_sub = zeros(n, 0);
    end

    partial_outputs{s} = out_sub;
    partial_P_joint{s} = P_sub;
    partial_V_joint{s} = V_sub;

    %% -----------------------------------------------------
    %% Diagnostics for target subsets 14 and 24.
    %% Important: use out_sub.blocks{local}.P_VAk, i.e. the
    %% initial low-rank row-subspace projector used inside this
    %% specific mssat call for the residualized subset.
    %% No second mssat call is made here.
    %% -----------------------------------------------------
    if target_mask(s)
        target_key = sprintf('subset_%s', lab);
        partial_projection_diag.(target_key) = struct();
        partial_projection_diag.(target_key).subset_label = lab;
        partial_projection_diag.(target_key).subset = idx;
        partial_projection_diag.(target_key).rank = r_hat_s;
        partial_projection_diag.(target_key).P_joint = P_sub;
        partial_projection_diag.(target_key).V_joint = V_sub;
        partial_projection_diag.(target_key).out_sub = out_sub;

        for a = 1:numel(idx)
            k = idx(a);
            %% block name without local helper function
            bname = sprintf('block%d', k);
            try
                if exist('fields', 'var') && ~isempty(fields)
                    if iscell(fields)
                        tmp_name = fields{k};
                    elseif isstring(fields)
                        tmp_name = fields(k);
                    else
                        tmp_name = [];
                    end

                    if iscell(tmp_name)
                        tmp_name = tmp_name{1};
                    end

                    if isstring(tmp_name)
                        bname = char(tmp_name);
                    elseif ischar(tmp_name)
                        bname = tmp_name;
                    end
                end
            catch
                bname = sprintf('block%d', k);
            end

            %% Initial low-rank projector for block k within this subset call.
            if isfield(out_sub, 'blocks') && numel(out_sub.blocks) >= a && ...
                    isfield(out_sub.blocks{a}, 'P_VAk') && ~isempty(out_sub.blocks{a}.P_VAk)
                P_init = out_sub.blocks{a}.P_VAk;
            else
                warning('P_VAk missing for subset %s, block %d. Using zero projector.', lab, k);
                P_init = zeros(n, n);
            end

            initial_rank = rank(P_init);

            if r_hat_s > 0 && ~isempty(V_sub)
                %% Subspace projection error:
                %% ||(I - P_init) V_J||_F^2 / ||V_J||_F^2.
                subspace_projection_error = norm((In - P_init) * V_sub, 'fro')^2 / ...
                    max(norm(V_sub, 'fro')^2, eps);

                %% Projected joint row subspace and principal angles.
                V_proj_raw = P_init * V_sub;
                [Qproj, Rproj] = qr(V_proj_raw, 0);

                if isempty(Rproj)
                    V_proj = zeros(n, 0);
                else
                    keep = abs(diag(Rproj)) > 1e-10;
                    V_proj = Qproj(:, keep);
                end

                if size(V_proj, 2) > 0
                    cosvals = svd(V_sub' * V_proj);
                    cosvals = min(max(cosvals, 0), 1);
                    theta_deg = acos(cosvals) * 180 / pi;
                else
                    theta_deg = NaN(r_hat_s, 1);
                end

                %% Pad if projection rank is deficient so rank-7 gives 7 entries.
                if numel(theta_deg) < r_hat_s
                    theta_deg = [theta_deg(:); NaN(r_hat_s - numel(theta_deg), 1)];
                else
                    theta_deg = theta_deg(:);
                end

                % %% Data-space matrix projection error.
                % %% Use the centered residual matrix, matching center=true inside mssat.
                % Xk_diag = bsxfun(@minus, R_list{k}, mean(R_list{k}, 2));
                % Jk = Xk_diag * P_sub;
                % Jk_proj = Jk * P_init;
                % 
                % matrix_projection_error = norm(Jk - Jk_proj, 'fro')^2 / ...
                %     max(norm(Jk, 'fro')^2, eps);

                %% Data-space matrix projection error.
                %% Use the current residual matrix directly, without additional centering.
                Xk_diag = R_list{k};

                Jk = Xk_diag * P_sub;
                Jk_proj = Jk * P_init;

                matrix_projection_error = norm(Jk - Jk_proj, 'fro')^2 / ...
                    max(norm(Jk, 'fro')^2, eps);
            else
                subspace_projection_error = NaN;
                matrix_projection_error = NaN;
                theta_deg = NaN;
                V_proj = zeros(n, 0);
                Jk = zeros(size(R_list{k}));
                Jk_proj = zeros(size(R_list{k}));
            end

            %% Save structured diagnostics.
            partial_projection_diag.(target_key).block(a).block_id = k;
            partial_projection_diag.(target_key).block(a).block_name = bname;
            partial_projection_diag.(target_key).block(a).initial_rank = initial_rank;
            partial_projection_diag.(target_key).block(a).P_initial = P_init;
            partial_projection_diag.(target_key).block(a).V_projected = V_proj;
            partial_projection_diag.(target_key).block(a).subspace_projection_error = subspace_projection_error;
            partial_projection_diag.(target_key).block(a).matrix_projection_error = matrix_projection_error;
            partial_projection_diag.(target_key).block(a).principal_angles_deg = theta_deg;
            partial_projection_diag.(target_key).block(a).J_matrix = Jk;
            partial_projection_diag.(target_key).block(a).J_matrix_projected = Jk_proj;

            %% Wide summary row.
            diag_subset{end+1, 1} = lab;
            diag_block_id(end+1, 1) = k;
            diag_block_name{end+1, 1} = bname;
            diag_joint_rank(end+1, 1) = r_hat_s;
            diag_initial_rank(end+1, 1) = initial_rank;
            diag_subspace_projection_error(end+1, 1) = subspace_projection_error;
            diag_matrix_projection_error(end+1, 1) = matrix_projection_error;
            diag_mean_angle_deg(end+1, 1) = mean(theta_deg, 'omitnan');
            diag_max_angle_deg(end+1, 1) = max(theta_deg, [], 'omitnan');
            diag_principal_angles_deg{end+1, 1} = theta_deg(:)';

            %% Long angle rows.
            for aa = 1:numel(theta_deg)
                angle_subset{end+1, 1} = lab;
                angle_block_id(end+1, 1) = k;
                angle_block_name{end+1, 1} = bname;
                angle_index(end+1, 1) = aa;
                angle_value_deg(end+1, 1) = theta_deg(aa);
            end
        end
    end

    %% Residualize current subset before moving to the next subset.
    if r_hat_s > 0 && ~isempty(P_sub)
        for kk = 1:numel(idx)
            j = idx(kk);
            R_list{j} = R_list{j} * (In - P_sub);
        end
    end
end

time.mssat_partial = toc;
fprintf('Total time (sequential partial joint): %.2f seconds\n', time.mssat_partial);

partial_ranks = cell2struct(num2cell(partial_ranks_vec), partial_rank_names, 2);
partial_rank_table = array2table(partial_ranks_vec, 'VariableNames', partial_rank_names);

disp('--- Partial joint ranks in specified order ---');
disp(partial_rank_table);

proj_tbl = table( ...
    diag_subset, ...
    diag_block_id, ...
    diag_block_name, ...
    diag_joint_rank, ...
    diag_initial_rank, ...
    diag_subspace_projection_error, ...
    diag_matrix_projection_error, ...
    diag_mean_angle_deg, ...
    diag_max_angle_deg, ...
    diag_principal_angles_deg, ...
    'VariableNames', { ...
        'subset', ...
        'block_id', ...
        'block_name', ...
        'joint_rank', ...
        'initial_rank', ...
        'subspace_projection_error', ...
        'matrix_projection_error', ...
        'mean_principal_angle_deg', ...
        'max_principal_angle_deg', ...
        'principal_angles_deg'});

angle_tbl = table( ...
    angle_subset, ...
    angle_block_id, ...
    angle_block_name, ...
    angle_index, ...
    angle_value_deg, ...
    'VariableNames', { ...
        'subset', ...
        'block_id', ...
        'block_name', ...
        'angle_index', ...
        'principal_angle_deg'});

disp('===== MSSAT partial projection diagnostic table =====');
disp(proj_tbl);

disp('===== MSSAT partial principal angle long table =====');
disp(angle_tbl);

%% =========================================================
%% Save MSSAT partial result and diagnostics
%% =========================================================
if exist(save_file, 'file') == 2
    S = load(save_file, 'save_struct');
    save_struct = S.save_struct;
else
    save_struct = struct();
end

save_struct.source_data_file = '\Matlab\AJIVE_Project\DataExample\TCGA.mat';
save_struct.realdata = realdata;
save_struct.fields = fields;
save_struct.time = time;
save_struct.jrank = jrank;

if exist('result_real', 'var')
    save_struct.result_real = result_real;
end

if ~isfield(save_struct, 'mssat') || ~isstruct(save_struct.mssat)
    save_struct.mssat = struct();
end

save_struct.mssat.was_run = true;
if exist('out_mssat', 'var')
    save_struct.mssat.out_mssat = out_mssat;
end
save_struct.mssat.partial_ranks = partial_ranks;
save_struct.mssat.partial_rank_table = partial_rank_table;
save_struct.mssat.partial_outputs = partial_outputs;
save_struct.mssat.partial_P_joint = partial_P_joint;
save_struct.mssat.partial_V_joint = partial_V_joint;
save_struct.mssat.partial_projection_diag = partial_projection_diag;
save_struct.mssat.partial_projection_table = proj_tbl;
save_struct.mssat.partial_projection_angle_table = angle_tbl;

if isfield(time, 'mssat')
    save_struct.mssat.time = time.mssat;
end
save_struct.mssat.time_partial = time.mssat_partial;

if isfield(jrank, 'mssat')
    save_struct.mssat.jrank = jrank.mssat;
end

if exist('finalTable', 'var')
    save_struct.finalTable = finalTable;
end

save(save_file, 'save_struct', '-v7.3');
fprintf('[Saved MSSAT partial outputs and projection diagnostics] %s\n', save_file);

%% =========================================================
%% Final message
%% =========================================================
fprintf('\n============================================\n');
fprintf('MSSAT partial projection diagnostics finished and saved.\n');
fprintf('File: %s\n', save_file);
fprintf('============================================\n');

end

%% =========================================================
%% Additional TCGA projection-error and principal-angle diagnostics
%% Computes:
%%   1) DIVAS 4-block global common structure vs initial row subspaces
%%   2) MSSAT CN-MUT pairwise common structure vs initial row subspaces
%% =========================================================
disp(' ');
disp('================ Projection-error and principal-angle diagnostics ================');

%% Reload the latest saved structure because the partial block above may have updated it.
if exist(save_file, 'file') == 2
    S = load(save_file, 'save_struct');
    save_struct = S.save_struct;
else
    error('save_file does not exist: %s', save_file);
end

if ~isfield(save_struct, 'realdata')
    error('save_struct.realdata is missing.');
end
if ~isfield(save_struct, 'fields')
    error('save_struct.fields is missing.');
end
if ~isfield(save_struct, 'divas') || ~isfield(save_struct.divas, 'DIVASout')
    error('save_struct.divas.DIVASout is missing. Run the DIVAS section first.');
end

realdata = save_struct.realdata;
fields = save_struct.fields;
DIVASout = save_struct.divas.DIVASout;

X_list = realdata;
K = numel(X_list);

if K ~= 4
    error('This diagnostics block assumes exactly 4 data blocks.');
end

field_str = string(fields);

relErr = @(A,B) (norm(A-B,'fro')^2) / max(norm(A,'fro')^2, eps);

%% ---------------------------------------------------------
%% Initial row subspaces: top right singular subspaces after row-centering.
%% ---------------------------------------------------------
Xc_list = cell(K,1);
V_init = cell(K,1);
P_init = cell(K,1);
r_init = zeros(K,1);

for k = 1:K
    Xk = double(X_list{k});
    Xc = bsxfun(@minus, Xk, mean(Xk, 2));
    Xc_list{k} = Xc;

    [~, rk] = BEMA_combined(Xc);
    r_init(k) = rk;

    [~,~,Vk] = svd(Xc, 'econ');
    V_init{k} = Vk(:, 1:rk);
    P_init{k} = V_init{k} * V_init{k}';
end

disp('Initial ranks from BEMA:');
disp(r_init');

%% =========================================================
%% PART A. DIVAS 4-block global common structure
%% =========================================================
disp(' ');
disp('=========================================================');
disp('PART A. DIVAS 4-block global common structure');
disp('=========================================================');

if ~isfield(DIVASout, 'keyIdxMap') || ~isfield(DIVASout, 'matBlocks')
    error('DIVASout does not contain keyIdxMap / matBlocks.');
end

keys_raw = DIVASout.keyIdxMap.keys;
vals_raw = DIVASout.keyIdxMap.values;

global_key = [];
global_pos = [];

for i = 1:numel(keys_raw)
    idxs = sort(vals_raw{i}(:))';
    if isequal(idxs, 1:K)
        global_key = keys_raw{i};
        global_pos = i;
        break;
    end
end

if isempty(global_pos)
    error('Could not find the 4-block global structure in DIVASout.keyIdxMap.');
end

fprintf('DIVAS global key = %s (position %d)\n', string(global_key), global_pos);

XJ_divas = cell(K,1);
XJ_divas_proj = cell(K,1);
divas_matrix_projection_error = zeros(K,1);
divas_common_rank = zeros(K,1);
divas_angles_deg = cell(K,1);
divas_largest_angle_deg = zeros(K,1);

for k = 1:K
    blk = DIVASout.matBlocks{k};

    if isstruct(blk) && isfield(blk, 'values') && isa(blk.values, 'containers.Map')
        try
            tmp = values(blk.values, {global_key});
            XJ_divas{k} = tmp{1};
        catch
            vals_k = blk.values;
            XJ_divas{k} = vals_k{global_pos};
        end
    elseif isstruct(blk) && isfield(blk, 'values') && iscell(blk.values)
        XJ_divas{k} = blk.values{global_pos};
    elseif isa(blk, 'containers.Map')
        tmp = values(blk, {global_key});
        XJ_divas{k} = tmp{1};
    else
        error('Unsupported DIVASout.matBlocks{%d} format.', k);
    end

    XJ_divas_proj{k} = XJ_divas{k} * P_init{k};
    divas_matrix_projection_error(k) = relErr(XJ_divas{k}, XJ_divas_proj{k});

    [~,~,VJk] = svd(XJ_divas{k}, 'econ');
    rJk = rank(XJ_divas{k});
    divas_common_rank(k) = rJk;

    if rJk > 0
        VJk = VJk(:, 1:rJk);
        cosvals = svd(V_init{k}' * VJk);
        cosvals = min(max(cosvals, 0), 1);
        divas_angles_deg{k} = acosd(cosvals);
        divas_largest_angle_deg(k) = subspace(V_init{k}, VJk) * 180 / pi;
    else
        divas_angles_deg{k} = [];
        divas_largest_angle_deg(k) = NaN;
    end
end

T_divas_projection_angle = table( ...
    (1:K)', ...
    field_str(:), ...
    r_init(:), ...
    divas_common_rank(:), ...
    divas_matrix_projection_error(:), ...
    divas_largest_angle_deg(:), ...
    divas_angles_deg(:), ...
    'VariableNames', { ...
        'Block', ...
        'Name', ...
        'InitialRank', ...
        'CommonRowRank', ...
        'MatrixProjectionError', ...
        'LargestSubspaceAngleDeg', ...
        'PrincipalAnglesDeg'} ...
);

disp('===== DIVAS projection-error and principal-angle summary =====');
disp(T_divas_projection_angle);

for k = 1:K
    fprintf('\n[DIVAS] Block %d (%s) principal angles (deg):\n', k, char(field_str(k)));
    disp(divas_angles_deg{k}');
end

%% =========================================================
%% PART B. MSSAT CN-MUT pairwise common structure
%% =========================================================
disp(' ');
disp('=========================================================');
disp('PART B. MSSAT CN-MUT pairwise common structure');
disp('=========================================================');

idx_CN = find(strcmpi(field_str, "CN"), 1);
idx_MUT = find(strcmpi(field_str, "MUT"), 1);

if isempty(idx_CN) || isempty(idx_MUT)
    fprintf('\n[Warning] Could not automatically find CN / MUT in fields.\n');
    disp(field_str);
    error('Please set idx_CN and idx_MUT manually.');
end

fprintf('CN block index  = %d\n', idx_CN);
fprintf('MUT block index = %d\n', idx_MUT);

X_pair = {Xc_list{idx_CN}, Xc_list{idx_MUT}};

if exist('mssat', 'file') == 2
    out_pair = mssat( ...
        X_pair, ...
        'center', false, ...
        'alpha', normcdf(-6), ...
        'rank_method', 'bema', ...
        'sigma_method', 'bema', ...
        'test_mode', 'alpha', ...
        'mode', 'test', ...
        'alt_angle_deg', 10, ...
        'verbose', false);
else
    error('mssat is not available on the MATLAB path.');
end

if ~isfield(out_pair, 'P_joint') || isempty(out_pair.P_joint)
    error('Pairwise MSSAT result does not contain P_joint.');
end

P_pair = out_pair.P_joint;

pair_block_ids = [idx_CN; idx_MUT];
pair_names = field_str(pair_block_ids);
pair_common_rank = zeros(2,1);
pair_matrix_projection_error = zeros(2,1);
pair_largest_angle_deg = zeros(2,1);
pair_angles_deg = cell(2,1);
XJ_pair = cell(2,1);
XJ_pair_proj = cell(2,1);

for a = 1:2
    k = pair_block_ids(a);

    XJ_pair{a} = Xc_list{k} * P_pair;
    XJ_pair_proj{a} = XJ_pair{a} * P_init{k};
    pair_matrix_projection_error(a) = relErr(XJ_pair{a}, XJ_pair_proj{a});

    [~,~,VJk] = svd(XJ_pair{a}, 'econ');
    rJk = rank(XJ_pair{a});
    pair_common_rank(a) = rJk;

    if rJk > 0
        VJk = VJk(:, 1:rJk);
        cosvals = svd(V_init{k}' * VJk);
        cosvals = min(max(cosvals, 0), 1);
        pair_angles_deg{a} = acosd(cosvals);
        pair_largest_angle_deg(a) = subspace(V_init{k}, VJk) * 180 / pi;
    else
        pair_angles_deg{a} = [];
        pair_largest_angle_deg(a) = NaN;
    end
end

T_mssat_pair_projection_angle = table( ...
    pair_block_ids(:), ...
    pair_names(:), ...
    r_init(pair_block_ids), ...
    pair_common_rank(:), ...
    pair_matrix_projection_error(:), ...
    pair_largest_angle_deg(:), ...
    pair_angles_deg(:), ...
    'VariableNames', { ...
        'Block', ...
        'Name', ...
        'InitialRank', ...
        'CommonRowRank', ...
        'MatrixProjectionError', ...
        'LargestSubspaceAngleDeg', ...
        'PrincipalAnglesDeg'} ...
);

disp('===== MSSAT CN-MUT projection-error and principal-angle summary =====');
disp(T_mssat_pair_projection_angle);

for a = 1:2
    fprintf('\n[MSSAT CN-MUT] Block %d (%s) principal angles (deg):\n', ...
        pair_block_ids(a), char(pair_names(a)));
    disp(pair_angles_deg{a}');
end

%% =========================================================
%% PART C. Combined comparison table
%% =========================================================
case_names = [ ...
    "DIVAS-" + field_str(:); ...
    "MSSAT-CN-MUT-" + pair_names(:) ...
];

matrix_projection_errors = [ ...
    divas_matrix_projection_error(:); ...
    pair_matrix_projection_error(:) ...
];

largest_angles = [ ...
    divas_largest_angle_deg(:); ...
    pair_largest_angle_deg(:) ...
];

T_projection_angle_compare = table( ...
    case_names, ...
    matrix_projection_errors, ...
    largest_angles, ...
    'VariableNames', { ...
        'Case', ...
        'MatrixProjectionError', ...
        'LargestSubspaceAngleDeg'} ...
);

disp('===== Combined projection-error and principal-angle comparison =====');
disp(T_projection_angle_compare);



