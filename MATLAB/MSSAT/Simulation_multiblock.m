%% =========================================================
%% Parallel 3-block simulation: AJIVE / DIVAS / MSSAT
%% Each worker uses one computational thread
%% Generated from simulation_multi.m
%% =========================================================
clear; clc;

run('~/MATLAB/cvx/cvx_setup.m');
run('~/MATLAB/cvx/cvx_startup.m');

%% ================== Paths ==================
addpath('~/MATLAB/MSSAT');
addpath('~/MATLAB/AJIVE_Project/AJIVECode');
addpath('~/MATLAB/Data-Integration-Via-Analysis-of-Subspaces/DJIVECode');
addpath('~/MATLAB/cvx');

%% ================== Parallel options ==================
ncores = 25;
maxNumCompThreads(1);

poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool('local', ncores);
elseif poolobj.NumWorkers ~= ncores
    delete(poolobj);
    poolobj = parpool('local', ncores);
end
pctRunOnAll maxNumCompThreads(1)

%% ================== User options ==================
% Available: 'ajive', 'divas', 'mssat'
methods = {'ajive', 'divas', 'mssat'};
rep = 100;
BASE_SEED = 20260514;

save_dir = '~/MATLAB/MSSAT/saved_results';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
raw_rep_dir = fullfile(save_dir, 'simulation_multi_parallel_raw');
if ~exist(raw_rep_dir, 'dir')
    mkdir(raw_rep_dir);
end

save_each_rep = true;
skip_existing_rep = true;  

%% ================== Simulation setting ==================
n = 500;
p1 = 400;
p2 = 300;
p3 = 200;

r1 = 6;
r2 = 9;
r3 = 12;

r = 3;
r12 = 1;
r13 = 1;
r23 = 2;

d = 1;
ang = pi/8;

d1 = d/2 * 3 * sqrt(n) * [6,5,4,7,3,2];
d2 = d/2 * 2 * sqrt(n) * [7,6,5,8,9,10,2,3,4];
d3 = d/2 * 4 * sqrt(n) * [10,9,8,7,5,6,2,3,4,11,12,13];

s1 = 1;
s2 = 1;
s3 = 1;

true_ranks = [r, r12, r13, r23, ...
              r1 - (r + r12 + r13), ...
              r2 - (r + r12 + r23), ...
              r3 - (r + r13 + r23)];

%% ================== Fixed bases shared across reps ==================
rng(0, 'twister');

v0 = random_orthonormal(n, r);
v_temp = random_orthonormal_orthogonal(r12 + r13 + r23, v0);
v12 = v_temp(:, 1:r12);
v13 = v_temp(:, (r12+1):(r12+r13));
v23 = v_temp(:, (r12+r13+1):(r12+r13+r23));

v11 = random_orthonormal_orthogonal(r1-r-r12-r13, [v0 v12 v13]);
v22 = random_orthonormal_orthogonal(r2-r-r12-r23, [v0 v12 v23]);
v33_temp = random_orthonormal_orthogonal(r3-r-r13-r23, [v0 v13 v23 v22]);
v33 = [v22*cos(ang)+v33_temp(:,1:(r2-r-r12-r23))*sin(ang), ...
       v33_temp(:,(r2-r-r12-r23+1):(r3-r-r13-r23))];

u1 = random_orthonormal(p1, r1);
u2 = random_orthonormal(p2, r2);
u3 = random_orthonormal(p3, r3);

fixed = struct();
fixed.n = n; fixed.p1 = p1; fixed.p2 = p2; fixed.p3 = p3;
fixed.r1 = r1; fixed.r2 = r2; fixed.r3 = r3;
fixed.r = r; fixed.r12 = r12; fixed.r13 = r13; fixed.r23 = r23;
fixed.d1 = d1; fixed.d2 = d2; fixed.d3 = d3;
fixed.s1 = s1; fixed.s2 = s2; fixed.s3 = s3;
fixed.ang = ang;
fixed.v0 = v0; fixed.v12 = v12; fixed.v13 = v13; fixed.v23 = v23;
fixed.v11 = v11; fixed.v22 = v22; fixed.v33 = v33;
fixed.u1 = u1; fixed.u2 = u2; fixed.u3 = u3;
fixed.true_ranks = true_ranks;
fixed.methods = methods;
fixed.BASE_SEED = BASE_SEED;

fprintf('\n================= Parallel simulation =================\n');
fprintf('methods = %s\n', strjoin(methods, ', '));
fprintf('rep = %d, ncores = %d\n', rep, ncores);
fprintf('raw_rep_dir = %s\n', raw_rep_dir);
fprintf('Each worker thread limit: maxNumCompThreads(1)\n');
fprintf('=======================================================\n\n');

%% ================== Preallocate result arrays ==================
method_template = struct();
for mm = 1:numel(methods)
    method = methods{mm};
    method_template.(method).rank = nan(rep, 7);
    method_template.(method).err  = nan(rep, 7);
    method_template.(method).time = nan(rep, 1);
    method_template.(method).error_msg = cell(rep, 1);
end

rep_results = cell(rep, 1);

overall_tic = tic;

%% ================== Parallel loop over rep ==================
parfor rep_idx = 1:rep
    maxNumCompThreads(1);

    out_file = fullfile(raw_rep_dir, sprintf('rep_%03d_result.mat', rep_idx));

    if skip_existing_rep && exist(out_file, 'file') == 2
        S = load(out_file, 'rep_result');
        rep_results{rep_idx} = S.rep_result;
        fprintf('[skip] rep %03d already exists\n', rep_idx);
        continue;
    end

    rep_result = run_one_rep_simulation_multi(rep_idx, fixed);

    if save_each_rep
        parsave_rep_result(out_file, rep_result);
    end

    rep_results{rep_idx} = rep_result;
    fprintf('[done] rep %03d / %03d, elapsed %.2fs\n', rep_idx, rep, rep_result.elapsed);
end

%% ================== Aggregate ==================
raw_results = method_template;

for rep_idx = 1:rep
    rr = rep_results{rep_idx};
    if isempty(rr)
        in_file = fullfile(raw_rep_dir, sprintf('rep_%03d_result.mat', rep_idx));
        if exist(in_file, 'file') == 2
            S = load(in_file, 'rep_result');
            rr = S.rep_result;
        else
            warning('Missing rep result: %d', rep_idx);
            continue;
        end
    end

    for mm = 1:numel(methods)
        method = methods{mm};
        raw_results.(method).rank(rep_idx, :) = rr.(method).rank;
        raw_results.(method).err(rep_idx, :)  = rr.(method).err;
        raw_results.(method).time(rep_idx)    = rr.(method).time;
        raw_results.(method).error_msg{rep_idx} = rr.(method).error_msg;
    end
end

raw_results.true_ranks = true_ranks;
raw_results.setting = fixed;
raw_results.total_time = toc(overall_tic);
raw_results.raw_rep_dir = raw_rep_dir;

%% ================== Summary tables ==================
structures = {'Global','J12','J13','J23','Ind1','Ind2','Ind3'}';
summary_map = struct();
time_map = struct();

for mm = 1:numel(methods)
    method = methods{mm};
    R = raw_results.(method).rank;
    E = raw_results.(method).err;

    mean_rank_hat = mean(R, 1, 'omitnan');
    prop_correct  = mean(R == true_ranks, 1, 'omitnan');
    mean_err      = mean(E, 1, 'omitnan');
    sd_err        = std(E, 0, 1, 'omitnan');

    summary_tbl = table( ...
        true_ranks(:), ...
        mean_rank_hat(:), ...
        prop_correct(:), ...
        mean_err(:), ...
        sd_err(:), ...
        'VariableNames', {'true_r','mean_r_hat','prop_correct','mean_relErr','sd_relErr'}, ...
        'RowNames', structures ...
    );

    summary_map.(method) = summary_tbl;
    time_map.(method) = sum(raw_results.(method).time, 'omitnan');

    fprintf('\n===== %s summary =====\n', upper(method));
    disp(summary_tbl);
    fprintf('%s total worker time sum: %.2f sec\n', upper(method), time_map.(method));
end

%% Rank summary table: mean rank (prop correct)
rank_rows = cell(numel(structures), 1 + numel(methods));
for ss = 1:numel(structures)
    rank_rows{ss,1} = structures{ss};
    for mm = 1:numel(methods)
        method = methods{mm};
        Tm = summary_map.(method);
        rank_rows{ss,mm+1} = sprintf('%.2f (%.2f)', Tm.mean_r_hat(ss), Tm.prop_correct(ss));
    end
end
rankTable = cell2table(rank_rows, 'VariableNames', ['Structure', methods]);

%% Error summary table: mean (sd)
err_rows = cell(numel(structures), 1 + numel(methods));
for ss = 1:numel(structures)
    err_rows{ss,1} = structures{ss};
    for mm = 1:numel(methods)
        method = methods{mm};
        Tm = summary_map.(method);
        err_rows{ss,mm+1} = sprintf('%.4f (%.4f)', Tm.mean_relErr(ss), Tm.sd_relErr(ss));
    end
end
errTable = cell2table(err_rows, 'VariableNames', ['Structure', methods]);

%% Time table
wall_time = raw_results.total_time;
time_vals = cell(1, numel(methods));
for mm = 1:numel(methods)
    method = methods{mm};
    time_vals{mm} = time_map.(method);
end
timeTable = cell2table(time_vals, 'VariableNames', methods, 'RowNames', {'sum_worker_time'});
timeTable.wall_time_sec = wall_time;

finalTable = struct();
finalTable.rank_multi3block = rankTable;
finalTable.err_multi3block = errTable;
finalTable.time_multi3block = timeTable;

fprintf('\n===== Combined rank table =====\n');
disp(rankTable);
fprintf('\n===== Combined error table =====\n');
disp(errTable);
fprintf('\n===== Time table =====\n');
disp(timeTable);

%% ================== Save final raw result ==================
save(fullfile(save_dir, 'sim_multi_result_parallel.mat'), ...
     'raw_results', 'summary_map', 'finalTable', '-v7.3');

fprintf('\nDONE. Saved:\n%s\n', fullfile(save_dir, 'sim_multi_result_parallel.mat'));

%% =========================================================
%% Local functions
%% =========================================================

function rep_result = run_one_rep_simulation_multi(rep_idx, fixed)
    maxNumCompThreads(1);

    methods = fixed.methods;
    n = fixed.n;
    p1 = fixed.p1; p2 = fixed.p2; p3 = fixed.p3;
    r1 = fixed.r1; r2 = fixed.r2; r3 = fixed.r3;
    r = fixed.r; r12 = fixed.r12; r13 = fixed.r13; r23 = fixed.r23;
    d1 = fixed.d1; d2 = fixed.d2; d3 = fixed.d3;
    s1 = fixed.s1; s2 = fixed.s2; s3 = fixed.s3;
    v0 = fixed.v0; v12 = fixed.v12; v13 = fixed.v13; v23 = fixed.v23;
    v11 = fixed.v11; v22 = fixed.v22; v33 = fixed.v33;
    u1 = fixed.u1; u2 = fixed.u2; u3 = fixed.u3;
    BASE_SEED = fixed.BASE_SEED;

    rep_result = struct();
    rep_result.rep_idx = rep_idx;
    rep_result.elapsed = NaN;

    for mm = 1:numel(methods)
        method = methods{mm};
        rep_result.(method).rank = nan(1,7);
        rep_result.(method).err = nan(1,7);
        rep_result.(method).time = NaN;
        rep_result.(method).error_msg = '';
    end

    rep_tic = tic;

    %% Generate one dataset per rep, shared by all methods for fair comparison
    rng(BASE_SEED + rep_idx, 'twister');

    rotG2  = random_orthonormal(r,r);
    rotG3  = random_orthonormal(r,r);
    rot12  = random_orthonormal(r12,r12);
    rot13  = random_orthonormal(r13,r13);
    rot23  = random_orthonormal(r23,r23);

    basis1 = [v0,         v12,        v13,         v11];
    basis2 = [v0*rotG2,   v12*rot12,  v23,         v22];
    basis3 = [v0*rotG3,   v13*rot13,  v23*rot23,   v33];

    S1 = u1 * diag(d1) * basis1';
    S2 = u2 * diag(d2) * basis2';
    S3 = u3 * diag(d3) * basis3';

    X1 = S1 + s1 * randn(p1, n);
    X2 = S2 + s2 * randn(p2, n);
    X3 = S3 + s3 * randn(p3, n);

    idx0_1  = 1:r;
    idx12_1 = r+1             : r+r12;
    idx13_1 = r+r12+1         : r+r12+r13;
    idxI1_1 = r+r12+r13+1     : r1;

    idx0_2  = 1:r;
    idx12_2 = r+1             : r+r12;
    idx23_2 = r+r12+1         : r+r12+r23;
    idxI2_2 = r+r12+r23+1     : r2;

    idx0_3  = 1:r;
    idx13_3 = r+1             : r+r13;
    idx23_3 = r+r13+1         : r+r13+r23;
    idxI3_3 = r+r13+r23+1     : r3;

    make_comp = @(u, d, basis, idx_keep) ...
        u * diag( d(:) .* ismember((1:length(d))', idx_keep) ) * basis';

    trueG1  = make_comp(u1, d1, basis1, idx0_1);
    trueG2  = make_comp(u2, d2, basis2, idx0_2);
    trueG3  = make_comp(u3, d3, basis3, idx0_3);

    true12_1 = make_comp(u1, d1, basis1, idx12_1);
    true12_2 = make_comp(u2, d2, basis2, idx12_2);

    true13_1 = make_comp(u1, d1, basis1, idx13_1);
    true13_3 = make_comp(u3, d3, basis3, idx13_3);

    true23_2 = make_comp(u2, d2, basis2, idx23_2);
    true23_3 = make_comp(u3, d3, basis3, idx23_3);

    trueI1 = make_comp(u1, d1, basis1, idxI1_1);
    trueI2 = make_comp(u2, d2, basis2, idxI2_2);
    trueI3 = make_comp(u3, d3, basis3, idxI3_3);

    relErr = @(A,B) (norm(A-B,'fro')^2) / max(norm(A,'fro')^2, eps);

    %% Run selected methods
    for mm = 1:numel(methods)
        method = methods{mm};
        t0 = tic;

        try
            switch lower(method)
                case 'ajive'
                    [rank_out, err_out] = run_ajive_hier3(X1, X2, X3, ...
                        trueG1, trueG2, trueG3, true12_1, true12_2, ...
                        true13_1, true13_3, true23_2, true23_3, ...
                        trueI1, trueI2, trueI3, relErr);

                case 'divas'
                    [rank_out, err_out] = run_divas_hier3(X1, X2, X3, ...
                        trueG1, trueG2, trueG3, true12_1, true12_2, ...
                        true13_1, true13_3, true23_2, true23_3, ...
                        trueI1, trueI2, trueI3, relErr);

                case 'mssat'
                    [rank_out, err_out] = run_mssat_hier3(X1, X2, X3, n, ...
                        trueG1, trueG2, trueG3, true12_1, true12_2, ...
                        true13_1, true13_3, true23_2, true23_3, ...
                        trueI1, trueI2, trueI3, relErr);

                otherwise
                    error('Unknown method: %s', method);
            end

            rep_result.(method).rank = rank_out;
            rep_result.(method).err = err_out;
            rep_result.(method).error_msg = '';

        catch ME
            rep_result.(method).rank = nan(1,7);
            rep_result.(method).err = nan(1,7);
            rep_result.(method).error_msg = ME.message;
        end

        rep_result.(method).time = toc(t0);
    end

    rep_result.elapsed = toc(rep_tic);
end

function [rank_out, err_out] = run_ajive_hier3(X1, X2, X3, trueG1, trueG2, trueG3, true12_1, true12_2, true13_1, true13_3, true23_2, true23_3, trueI1, trueI2, trueI3, relErr)
    n = size(X1,2);
    datablock = {X1, X2, X3};
    dataname  = {'X1','X2','X3'};

    [~, r1_hat0] = BEMA_combined(X1);
    [~, r2_hat0] = BEMA_combined(X2);
    [~, r3_hat0] = BEMA_combined(X3);
    vecr0 = [r1_hat0, r2_hat0, r3_hat0];

    param0 = struct('dataname', {dataname}, ...
                    'iplot',   [0 0], ...
                    'ioutput', [0 0 0 0 0 0 1 1 0]);

    out0 = AJIVEMainMJ(datablock, vecr0, param0);

    G1_hat = out0.MatrixJoint{1};
    G2_hat = out0.MatrixJoint{2};
    G3_hat = out0.MatrixJoint{3};
    rG_hat = out0.rjoint;

    R1 = X1 - G1_hat;
    R2 = X2 - G2_hat;
    R3 = X3 - G3_hat;

    db12 = {R1, R2};
    [~, r1_12_hat0] = BEMA_combined(R1);
    [~, r2_12_hat0] = BEMA_combined(R2);
    param12 = struct('dataname', {{'R1','R2'}}, 'iplot', [0 0], 'ioutput', [0 0 0 0 0 0 1 1 0]);
    out12 = AJIVEMainMJ(db12, [r1_12_hat0, r2_12_hat0], param12);
    J12_1_hat = out12.MatrixJoint{1};
    J12_2_hat = out12.MatrixJoint{2};
    r12_hat = out12.rjoint;
    R1 = R1 - J12_1_hat;
    R2 = R2 - J12_2_hat;

    db13 = {R1, R3};
    [~, r1_13_hat0] = BEMA_combined(R1);
    [~, r3_13_hat0] = BEMA_combined(R3);
    param13 = struct('dataname', {{'R1','R3'}}, 'iplot', [0 0], 'ioutput', [0 0 0 0 0 0 1 1 0]);
    out13 = AJIVEMainMJ(db13, [r1_13_hat0, r3_13_hat0], param13);
    J13_1_hat = out13.MatrixJoint{1};
    J13_3_hat = out13.MatrixJoint{2};
    r13_hat = out13.rjoint;
    R1 = R1 - J13_1_hat;
    R3 = R3 - J13_3_hat;

    db23 = {R2, R3};
    [~, r2_23_hat0] = BEMA_combined(R2);
    [~, r3_23_hat0] = BEMA_combined(R3);
    param23 = struct('dataname', {{'R2','R3'}}, 'iplot', [0 0], 'ioutput', [0 0 0 0 0 0 1 1 0]);
    out23 = AJIVEMainMJ(db23, [r2_23_hat0, r3_23_hat0], param23);
    J23_2_hat = out23.MatrixJoint{1};
    J23_3_hat = out23.MatrixJoint{2};
    r23_hat = out23.rjoint;
    R2 = R2 - J23_2_hat;
    R3 = R3 - J23_3_hat;

    rI1_hat = max(r1_hat0 - rG_hat - r12_hat - r13_hat, 0);
    rI2_hat = max(r2_hat0 - rG_hat - r12_hat - r23_hat, 0);
    rI3_hat = max(r3_hat0 - rG_hat - r13_hat - r23_hat, 0);

    I1_hat = local_lowrank(R1, rI1_hat);
    I2_hat = local_lowrank(R2, rI2_hat);
    I3_hat = local_lowrank(R3, rI3_hat);

    eG  = mean([relErr(trueG1, G1_hat), relErr(trueG2, G2_hat), relErr(trueG3, G3_hat)]);
    e12 = mean([relErr(true12_1, J12_1_hat), relErr(true12_2, J12_2_hat)]);
    e13 = mean([relErr(true13_1, J13_1_hat), relErr(true13_3, J13_3_hat)]);
    e23 = mean([relErr(true23_2, J23_2_hat), relErr(true23_3, J23_3_hat)]);
    eI1 = relErr(trueI1, I1_hat);
    eI2 = relErr(trueI2, I2_hat);
    eI3 = relErr(trueI3, I3_hat);

    rank_out = [rG_hat, r12_hat, r13_hat, r23_hat, rI1_hat, rI2_hat, rI3_hat];
    err_out  = [eG, e12, e13, e23, eI1, eI2, eI3];
end

function [rank_out, err_out] = run_divas_hier3(X1, X2, X3, trueG1, trueG2, trueG3, true12_1, true12_2, true13_1, true13_3, true23_2, true23_3, trueI1, trueI2, trueI3, relErr)
    p1 = size(X1,1); p2 = size(X2,1); p3 = size(X3,1); n = size(X1,2);
    data_sim = {X1, X2, X3};
    paramstruct_divas = struct();
    DIVASout = DJIVEMainJP(data_sim, paramstruct_divas);

    G_hat  = {zeros(p1,n), zeros(p2,n), zeros(p3,n)};
    J12_hat = {zeros(p1,n), zeros(p2,n)};
    J13_hat = {zeros(p1,n), zeros(p3,n)};
    J23_hat = {zeros(p2,n), zeros(p3,n)};
    I_hat  = {zeros(p1,n), zeros(p2,n), zeros(p3,n)};

    keys_cell = DIVASout.keyIdxMap.keys;
    for kk = 1:numel(keys_cell)
        key_str = keys_cell{kk};
        blk_idx = DIVASout.keyIdxMap(key_str);
        blk_idx = blk_idx(:)';

        for bb = blk_idx
            M = local_get_divas_matrix(DIVASout.matBlocks{bb}, key_str);

            if numel(blk_idx) == 3
                G_hat{bb} = G_hat{bb} + M;
            elseif numel(blk_idx) == 2
                sblk = sort(blk_idx);
                if isequal(sblk, [1 2])
                    if bb == 1, J12_hat{1} = J12_hat{1} + M; else, J12_hat{2} = J12_hat{2} + M; end
                elseif isequal(sblk, [1 3])
                    if bb == 1, J13_hat{1} = J13_hat{1} + M; else, J13_hat{2} = J13_hat{2} + M; end
                elseif isequal(sblk, [2 3])
                    if bb == 2, J23_hat{1} = J23_hat{1} + M; else, J23_hat{2} = J23_hat{2} + M; end
                end
            elseif numel(blk_idx) == 1
                I_hat{bb} = I_hat{bb} + M;
            end
        end
    end

    rG_hat  = rank([G_hat{1}; G_hat{2}; G_hat{3}]);
    r12_hat = rank([J12_hat{1}; J12_hat{2}]);
    r13_hat = rank([J13_hat{1}; J13_hat{2}]);
    r23_hat = rank([J23_hat{1}; J23_hat{2}]);
    rI1_hat = rank(I_hat{1});
    rI2_hat = rank(I_hat{2});
    rI3_hat = rank(I_hat{3});

    eG  = mean([relErr(trueG1, G_hat{1}), relErr(trueG2, G_hat{2}), relErr(trueG3, G_hat{3})]);
    e12 = mean([relErr(true12_1, J12_hat{1}), relErr(true12_2, J12_hat{2})]);
    e13 = mean([relErr(true13_1, J13_hat{1}), relErr(true13_3, J13_hat{2})]);
    e23 = mean([relErr(true23_2, J23_hat{1}), relErr(true23_3, J23_hat{2})]);
    eI1 = relErr(trueI1, I_hat{1});
    eI2 = relErr(trueI2, I_hat{2});
    eI3 = relErr(trueI3, I_hat{3});

    rank_out = [rG_hat, r12_hat, r13_hat, r23_hat, rI1_hat, rI2_hat, rI3_hat];
    err_out  = [eG, e12, e13, e23, eI1, eI2, eI3];
end

function [rank_out, err_out] = run_mssat_hier3(X1, X2, X3, n, trueG1, trueG2, trueG3, true12_1, true12_2, true13_1, true13_3, true23_2, true23_3, trueI1, trueI2, trueI3, relErr)
    [s1_hat2, r1_hat] = BEMA_combined(X1);
    [s2_hat2, r2_hat] = BEMA_combined(X2);
    [s3_hat2, r3_hat] = BEMA_combined(X3);

    outG = mssat({X1, X2, X3}, ...
        'center', false, ...
        'alpha', normcdf(-6), ...
        'rank_method', 'given', ...
        'sigma_method', 'given', ...
        'r_init', [r1_hat, r2_hat, r3_hat], ...
        'sigma_init', [s1_hat2, s2_hat2, s3_hat2], ...
        'test_mode', 'alpha', ...
        'mode', 'test', ...
        'alt_angle_deg', 5, ...
        'verbose', false);

    rG_hat = double(outG.joint_rank);
    PG = outG.P_joint;
    I_n = eye(n);

    G1_hat = X1 * PG;
    G2_hat = X2 * PG;
    G3_hat = X3 * PG;

    R1 = X1 * (I_n - PG);
    R2 = X2 * (I_n - PG);
    R3 = X3 * (I_n - PG);

    [s1_12, r1_12] = BEMA_combined(R1);
    [s2_12, r2_12] = BEMA_combined(R2);
    out12 = mssat({R1, R2}, ...
        'center', false, 'alpha', normcdf(-6), ...
        'rank_method','given', 'sigma_method','given', ...
        'r_init', [r1_12, r2_12], 'sigma_init',[s1_12, s2_12], ...
        'test_mode', 'alpha', 'mode', 'test', 'alt_angle_deg', 5, 'verbose', false);
    r12_hat = double(out12.joint_rank);
    P12 = out12.P_joint;
    J12_1_hat = R1 * P12;
    J12_2_hat = R2 * P12;
    R1 = R1 * (I_n - P12);
    R2 = R2 * (I_n - P12);

    [s1_13, r1_13] = BEMA_combined(R1);
    [s3_13, r3_13] = BEMA_combined(R3);
    out13 = mssat({R1, R3}, ...
        'center', false, 'alpha', normcdf(-6), ...
        'rank_method','given', 'sigma_method','given', ...
        'r_init', [r1_13, r3_13], 'sigma_init',[s1_13, s3_13], ...
        'test_mode', 'alpha', 'mode', 'test', 'alt_angle_deg', 5, 'verbose', false);
    r13_hat = double(out13.joint_rank);
    P13 = out13.P_joint;
    J13_1_hat = R1 * P13;
    J13_3_hat = R3 * P13;
    R1 = R1 * (I_n - P13);
    R3 = R3 * (I_n - P13);

    [s2_23, r2_23] = BEMA_combined(R2);
    [s3_23, r3_23] = BEMA_combined(R3);
    out23 = mssat({R2, R3}, ...
        'center', false, 'alpha', normcdf(-6), ...
        'rank_method','given', 'sigma_method','given', ...
        'r_init', [r2_23, r3_23], 'sigma_init',[s2_23, s3_23], ...
        'test_mode', 'alpha', 'mode', 'test', 'alt_angle_deg', 5, 'verbose', false);
    r23_hat = double(out23.joint_rank);
    P23 = out23.P_joint;
    J23_2_hat = R2 * P23;
    J23_3_hat = R3 * P23;
    R2 = R2 * (I_n - P23);
    R3 = R3 * (I_n - P23);

    rI1_hat = r1_hat - rG_hat - r12_hat - r13_hat;
    rI2_hat = r2_hat - rG_hat - r12_hat - r23_hat;
    rI3_hat = r3_hat - rG_hat - r13_hat - r23_hat;

    I1_hat = local_lowrank(R1, rI1_hat);
    I2_hat = local_lowrank(R2, rI2_hat);
    I3_hat = local_lowrank(R3, rI3_hat);

    eG  = mean([relErr(trueG1, G1_hat), relErr(trueG2, G2_hat), relErr(trueG3, G3_hat)]);
    e12 = mean([relErr(true12_1, J12_1_hat), relErr(true12_2, J12_2_hat)]);
    e13 = mean([relErr(true13_1, J13_1_hat), relErr(true13_3, J13_3_hat)]);
    e23 = mean([relErr(true23_2, J23_2_hat), relErr(true23_3, J23_3_hat)]);
    eI1 = relErr(trueI1, I1_hat);
    eI2 = relErr(trueI2, I2_hat);
    eI3 = relErr(trueI3, I3_hat);

    rank_out = [rG_hat, r12_hat, r13_hat, r23_hat, rI1_hat, rI2_hat, rI3_hat];
    err_out  = [eG, e12, e13, e23, eI1, eI2, eI3];
end

function M = local_get_divas_matrix(Mk, key_str)
    if isa(Mk, 'containers.Map')
        M = Mk(key_str);
    elseif isstruct(Mk) && isfield(Mk, 'values')
        vals = Mk.values;
        key_num = str2double(key_str);
        if iscell(vals)
            M = vals{key_num};
        else
            M = vals(key_num);
        end
    else
        error('Unknown DIVAS matBlocks structure.');
    end
end

function Xhat = local_lowrank(X, rr)
    rr = max(round(rr), 0);
    if rr > 0
        rr = min(rr, min(size(X)));
        [U,S,V] = svd(X, 'econ');
        Xhat = U(:,1:rr) * S(1:rr,1:rr) * V(:,1:rr)';
    else
        Xhat = zeros(size(X));
    end
end

function parsave_rep_result(filename, rep_result)
    save(filename, 'rep_result', '-v7');
end


%% =========================================================
%% Combined result table with reconstruction error
%% =========================================================

methods = {'ajive', 'divas', 'mssat'};

structures = {'Global','J12','J13','J23','Ind1','Ind2','Ind3'}';

% ---------------------------------------------------------
% 0. summary table
% ---------------------------------------------------------
summary_map = struct();
time_map    = struct();

if exist('summary_tbl_ajive', 'var')
    summary_map.ajive = summary_tbl_ajive;
end
if exist('summary_tbl_divas', 'var')
    summary_map.divas = summary_tbl_divas;
end
if exist('summary_tbl', 'var')
    summary_map.mssat = summary_tbl;
end

if exist('time_sim', 'var') && isfield(time_sim, 'ajive_3block')
    time_map.ajive = time_sim.ajive_3block;
end
if exist('time_sim', 'var') && isfield(time_sim, 'divas_3block')
    time_map.divas = time_sim.divas_3block;
end
if exist('time_sim', 'var') && isfield(time_sim, 'mssat_3block')
    time_map.mssat = time_sim.mssat_3block;
end

% ---------------------------------------------------------
methods_exist = {};
for m = 1:numel(methods)
    if isfield(summary_map, methods{m})
        methods_exist{end+1} = methods{m}; 
    else
        fprintf('[skip] %s: summary table does not exist.\n', methods{m});
    end
end

if isempty(methods_exist)
    error('No available methods found.');
end

fprintf('Methods used in combined table: %s\n', strjoin(methods_exist, ', '));

% ---------------------------------------------------------
rank_rows = cell(numel(structures), 1 + numel(methods_exist));

for s = 1:numel(structures)
    rank_rows{s,1} = structures{s};

    for m = 1:numel(methods_exist)
        method = methods_exist{m};
        Tm = summary_map.(method);

        rank_rows{s,m+1} = sprintf('%.2f (%.2f)', ...
            Tm.mean_r_hat(s), Tm.prop_correct(s));
    end
end

rankTable = cell2table(rank_rows, ...
    'VariableNames', ['Structure', methods_exist]);

disp('===== Rank summary: mean rank (prop. correct) =====');
disp(rankTable);

% ---------------------------------------------------------
err_rows = cell(numel(structures), 1 + numel(methods_exist));

for s = 1:numel(structures)
    err_rows{s,1} = structures{s};

    for m = 1:numel(methods_exist)
        method = methods_exist{m};
        Tm = summary_map.(method);

        err_rows{s,m+1} = sprintf('%.4f (%.4f)', ...
            Tm.mean_relErr(s), Tm.sd_relErr(s));
    end
end

errTable = cell2table(err_rows, ...
    'VariableNames', ['Structure', methods_exist]);

disp('===== Reconstruction error summary: mean (sd) =====');
disp(errTable);

% ---------------------------------------------------------
time_vals = cell(1, numel(methods_exist));
for m = 1:numel(methods_exist)
    method = methods_exist{m};
    if isfield(time_map, method)
        time_vals{m} = time_map.(method);
    else
        time_vals{m} = NaN;
    end
end

timeTable = cell2table(time_vals, ...
    'VariableNames', methods_exist, ...
    'RowNames', {'time'});

disp('===== Time summary =====');
disp(timeTable);

% ---------------------------------------------------------
finalTable.rank_multi3block = rankTable;
finalTable.err_multi3block  = errTable;
finalTable.time_multi3block = timeTable;

disp('===== Combined simulation result table saved to finalTable =====');