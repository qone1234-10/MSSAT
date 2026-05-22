% AJIVE / DIVAS / MSSAT package on TCGA 4-block data

var_keep = {'result_real', 'finalTable'};
vars = setdiff(who, var_keep);
clear(vars{:});
clc

addpath('~\Matlab\MSSAT');
addpath('~\Matlab\AJIVE_Project\AJIVECode');
addpath('~\Matlab\Data-Integration-Via-Analysis-of-Subspaces\DJIVECode');

%% =========================================================
%% Load TCGA data
%% =========================================================
data_file = '~\Matlab\AJIVE_Project\DataExample\TCGA.mat';
datatmp = load(data_file);

realdata = datatmp.datablock;
fields   = datatmp.dataname;
X_list   = realdata(:).';
K        = numel(X_list);

%% =========================================================
%% Save setup
%% =========================================================
save_dir = '~\Matlab\MSSAT\saved_results';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

save_file = fullfile(save_dir, 'tcga_4block_analysis.mat');
source_data_file = data_file;

fprintf('Save file = %s\n', save_file);

%% initialize summary containers
time = struct();
jrank = struct();

%% =========================================================
%% AJIVE
%% =========================================================
disp('================ AJIVE ================');

datablock = X_list;
dataname = {'X', 'Y', 'Z', 'W'};

disp('Running AJIVE on TCGA 4-block data...');

tic

[~, r1_hat] = BEMA_combined(X_list{1});
[~, r2_hat] = BEMA_combined(X_list{2});
[~, r3_hat] = BEMA_combined(X_list{3});
[~, r4_hat] = BEMA_combined(X_list{4});

vecr = [r1_hat, r2_hat, r3_hat, r4_hat];
disp('Initial ranks from BEMA:');
disp(vecr);

paramstruct0 = struct('dataname', {dataname}, ...
                      'iplot', [0 0]);

outstruct0 = AJIVEMainMJ(datablock, vecr, paramstruct0);

time.ajive = toc;
jrank.ajive = [ ...
    outstruct0.rjoint, ...
    r1_hat - outstruct0.rjoint, ...
    r2_hat - outstruct0.rjoint, ...
    r3_hat - outstruct0.rjoint, ...
    r4_hat - outstruct0.rjoint ...
];

disp(time);

% ----- save AJIVE result immediately -----
save_struct = load_save_struct(save_file);
save_struct = attach_common_fields(save_struct, source_data_file, realdata, fields, time, jrank);

save_struct.ajive = struct();
save_struct.ajive.was_run = true;
save_struct.ajive.outstruct0 = outstruct0;
save_struct.ajive.initial_ranks = vecr;
save_struct.ajive.time = time.ajive;
save_struct.ajive.jrank = jrank.ajive;

if exist('finalTable', 'var')
    save_struct.finalTable = finalTable;
end

save(save_file, 'save_struct', '-v7.3');
fprintf('[Saved after AJIVE] %s\n', save_file);

%% =========================================================
%% DIVAS
%% =========================================================
disp('================ DIVAS ================');

tic

DIVASout = DJIVEMainJP(X_list);

time.divas = toc;
jrank.divas = [ ...
    DIVASout.rjoint{1}(1), ...
    DIVASout.rjoint{1}(2), ...
    DIVASout.rjoint{2}(2), ...
    DIVASout.rjoint{3}(2), ...
    DIVASout.rjoint{4}(2) ...
];

disp(time);

% ----- save DIVAS result immediately -----
save_struct = load_save_struct(save_file);
save_struct = attach_common_fields(save_struct, source_data_file, realdata, fields, time, jrank);

save_struct.divas = struct();
save_struct.divas.was_run = true;
save_struct.divas.DIVASout = DIVASout;
save_struct.divas.time = time.divas;
save_struct.divas.jrank = jrank.divas;

if exist('finalTable', 'var')
    save_struct.finalTable = finalTable;
end

save(save_file, 'save_struct', '-v7.3');
fprintf('[Saved after DIVAS] %s\n', save_file);


%% =========================================================
%% MSSAT
%% =========================================================
disp('================ MSSAT ================');

fprintf('Running MSSAT on %d data blocks (BRCA)...\n', K);

tic
out_mssat = mssat(X_list, ...
    'center', true, ...
    'alpha', normcdf(-6), ...
    'rank_method', 'bema', ...
    'sigma_method', 'bema', ...
    'test_mode', 'alpha', ...
    'alt_angle_deg', 10, ...
    'mode', 'test', ...
    'verbose', true);
time.mssat = toc;

% ======================= stable rank parsing =======================
if isfield(out_mssat, 'blocks') && ~isempty(out_mssat.blocks)
    r_individual = zeros(1, K);
    for k = 1:K
        try
            Vj_block = out_mssat.blocks{k}.V_joint_block;
            Vind     = out_mssat.blocks{k}.V_indiv;
            r_individual(k) = size(Vj_block, 2) + size(Vind, 2);
        catch
            try
                Vk_block = out_mssat.blocks{k}.P_VAk;
                r_individual(k) = rank(Vk_block);
            catch
                r_individual(k) = 0;
            end
        end
    end
elseif isfield(out_mssat, 'options') && ...
       isfield(out_mssat.options, 'r_init') && ...
       ~isempty(out_mssat.options.r_init)
    r_individual = double(out_mssat.options.r_init);
else
    r_individual = zeros(1, K);
    for k = 1:K
        [~, rk_tmp] = BEMA_combined(X_list{k});
        r_individual(k) = rk_tmp;
    end
end

rJ = double(out_mssat.joint_rank);
r_individual = double(r_individual);
r_specific = max(r_individual - rJ, 0);
r_specific = double(r_specific);

jrank.mssat = [rJ, r_specific];

disp('-------------------------------------------------------');
fprintf('>>> MSSAT finished.\n');
fprintf('Estimated joint rank: %d\n', rJ);
fprintf('Individual ranks (after joint removal): [%s]\n', num2str(r_specific));
disp('-------------------------------------------------------');

% ----- save MSSAT result immediately -----
save_struct = load_save_struct(save_file);
save_struct = attach_common_fields(save_struct, source_data_file, realdata, fields, time, jrank);

save_struct.mssat = struct();
save_struct.mssat.was_run = true;
save_struct.mssat.out_mssat = out_mssat;
save_struct.mssat.time = time.mssat;
save_struct.mssat.jrank = jrank.mssat;

if exist('finalTable', 'var')
    save_struct.finalTable = finalTable;
end

save(save_file, 'save_struct', '-v7.3');
fprintf('[Saved after MSSAT] %s\n', save_file);

%% =========================================================
%% Summary table
%% =========================================================
result_real.tcga = combineRankAndTime(jrank, time);

disp("result from tcga");
disp(result_real.tcga);

% ----- save summary immediately -----
save_struct = load_save_struct(save_file);
save_struct = attach_common_fields(save_struct, source_data_file, realdata, fields, time, jrank);
save_struct.result_real = result_real;

if exist('finalTable', 'var')
    save_struct.finalTable = finalTable;
end

save(save_file, 'save_struct', '-v7.3');
fprintf('[Saved after summary] %s\n', save_file);

%% =========================================================
%% MSSAT partially joint
%% =========================================================
disp('================ MSSAT partial ================');

fprintf('Running sequential partially joint MSSAT on %d data blocks (BRCA)...\n', K);

% row-centered residual initialization
R_list = cell(1, K);
for k = 1:K
    Xk = X_list{k};
    R_list{k} = bsxfun(@minus, Xk, mean(Xk, 2));
end

n  = size(R_list{1}, 2);
In = eye(n);

subset_list = { [1 2 3], [1 2 4], [1 3 4], [2 3 4], ...
                [1 2],   [1 3],   [1 4],   [2 3],   [2 4],   [3 4] };
subset_labels = {'123','124','134','234','12','13','14','23','24','34'};

partial_ranks_vec = zeros(1, numel(subset_list));

fprintf('\n===== Sequential partially joint rank estimation (4 blocks, BRCA) =====\n');
tic;
for s = 1:numel(subset_list)
    idx = subset_list{s};
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

    fprintf('Subset %s: estimated joint rank = %d\n', subset_labels{s}, r_hat_s);

    P_sub = out_sub.P_joint;
    if r_hat_s > 0 && ~isempty(P_sub)
        for kk = 1:numel(idx)
            j = idx(kk);
            R_list{j} = R_list{j} * (In - P_sub);
        end
    end
end
time.mssat_partial = toc;
fprintf('Total time (sequential partial joint): %.2f seconds\n', time.mssat_partial);

partial_rank_names = {'r123','r124','r134','r234','r12','r13','r14','r23','r24','r34'};
partial_ranks = cell2struct(num2cell(partial_ranks_vec), partial_rank_names, 2);

partial_rank_table = array2table(partial_ranks_vec, ...
    'VariableNames', partial_rank_names);

disp('--- Partial joint ranks in specified order ---');
disp(partial_rank_table);

% ----- save MSSAT partial result immediately -----
save_struct = load_save_struct(save_file);
save_struct = attach_common_fields(save_struct, source_data_file, realdata, fields, time, jrank);
save_struct.result_real = result_real;

if ~isfield(save_struct, 'mssat') || ~isstruct(save_struct.mssat)
    save_struct.mssat = struct();
end

save_struct.mssat.was_run = true;
if exist('out_mssat', 'var')
    save_struct.mssat.out_mssat = out_mssat;
end
save_struct.mssat.partial_ranks = partial_ranks;
save_struct.mssat.partial_rank_table = partial_rank_table;
save_struct.mssat.time = time.mssat;
save_struct.mssat.time_partial = time.mssat_partial;
save_struct.mssat.jrank = jrank.mssat;

if exist('finalTable', 'var')
    save_struct.finalTable = finalTable;
end

save(save_file, 'save_struct', '-v7.3');
fprintf('[Saved after MSSAT partial] %s\n', save_file);

%% =========================================================
%% Final message
%% =========================================================
fprintf('\n============================================\n');
fprintf('TCGA 4-block analysis finished and saved.\n');
fprintf('File: %s\n', save_file);
fprintf('============================================\n');

%% =========================================================
%% Local helper functions
%% =========================================================
function save_struct = load_save_struct(save_file)
    if exist(save_file, 'file')
        S = load(save_file, 'save_struct');
        if isfield(S, 'save_struct')
            save_struct = S.save_struct;
        else
            save_struct = struct();
        end
    else
        save_struct = struct();
    end
end

function save_struct = attach_common_fields(save_struct, source_data_file, realdata, fields, time, jrank)
    save_struct.source_data_file = source_data_file;
    save_struct.realdata = realdata;
    save_struct.fields = fields;
    save_struct.time = time;
    save_struct.jrank = jrank;
end
