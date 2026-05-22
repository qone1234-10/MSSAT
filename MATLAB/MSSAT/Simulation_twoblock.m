%% =========================================================
%% Parallel simulation: 2-block grid, 25 cores
%% AJIVE / DIVAS / MSSAT
%% rep-level parallelization
%% =========================================================

clear; clc;

run('~/MATLAB/cvx/cvx_setup.m');
run('~/MATLAB/cvx/cvx_startup.m');

addpath('~/MATLAB/MSSAT');
addpath('~/MATLAB/AJIVE_Project/AJIVECode');
addpath('~/MATLAB/Data-Integration-Via-Analysis-of-Subspaces/DJIVECode');
addpath('~/MATLAB/cvx');
%% options
methods = {'ajive', 'divas', 'mssat'};

rep = 100;
ncores = 50;

BASE_SEED = 2026;

raw_full_dir = '~/MATLAB/sim_dang_raw_full_parallel';
if ~exist(raw_full_dir, 'dir')
    mkdir(raw_full_dir);
end

%% basic setting
n  = 500;
p1 = 400;
p2 = 300;

r1 = 6;
r2 = 9;
r  = 3;

s1 = 1;
s2 = 1;

true_r = r;
r1_ind_true = r1 - true_r;
r2_ind_true = r2 - true_r;

%% grid
d_vals   = [0.1 0.25 0.5 1 2];
ang_vals = [0 pi/8 pi/4 3*pi/8];

nd = numel(d_vals);
na = numel(ang_vals);
nm = numel(methods);

%% fixed left singular bases
rng(0);
u1 = random_orthonormal(p1, r1);
u2 = random_orthonormal(p2, r2);

fprintf('\nParallel simulation start\n');
fprintf('rep = %d, cores = %d\n', rep, ncores);
fprintf('methods = %s\n', strjoin(methods, ', '));
fprintf('raw dir = %s\n\n', raw_full_dir);

%% parallel pool
maxNumCompThreads(1);

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local', ncores);
elseif poolobj.NumWorkers ~= ncores
    delete(poolobj);
    parpool('local', ncores);
end

%% =========================================================
%% parfor over rep
%% =========================================================

parfor rep_idx = 1:rep

    %% worker path setup

    maxNumCompThreads(1);
    addpath('~/MATLAB/MSSAT');
    addpath('~/MATLAB/AJIVE_Project/AJIVECode');
    addpath('~/MATLAB/Data-Integration-Via-Analysis-of-Subspaces/DJIVECode');
    addpath('~/MATLAB/cvx');

    rep_tic = tic;

    rep_result = struct();
    rep_result.rep_idx = rep_idx;
    rep_result.d_vals = d_vals;
    rep_result.ang_vals = ang_vals;
    rep_result.methods = methods;

    % Skip this replication only if all requested method-specific files already exist.
    all_method_files_exist = true;
    for mm_check = 1:nm
        method_check = methods{mm_check};
        method_file_check = fullfile(raw_full_dir, ...
            sprintf('rep_%03d_%s_result.mat', rep_idx, method_check));
        if exist(method_file_check, 'file') ~= 2
            all_method_files_exist = false;
            break;
        end
    end

    if all_method_files_exist
        fprintf('[skip] rep %03d: all requested method files already exist\n', rep_idx);
        continue;
    end

    for mm = 1:nm
        method = methods{mm};

        rep_result.(method).rJ = nan(nd, na);
        rep_result.(method).rI1 = nan(nd, na);
        rep_result.(method).rI2 = nan(nd, na);

        rep_result.(method).init_rank1 = nan(nd, na);
        rep_result.(method).init_rank2 = nan(nd, na);

        rep_result.(method).time = nan(nd, na);

        rep_result.(method).err_overall = nan(nd, na, 2);
        rep_result.(method).err_joint = nan(nd, na, 2);
        rep_result.(method).err_ind = nan(nd, na, 2);

        rep_result.(method).error_msg = cell(nd, na);
    end

    for id = 1:nd
        d = d_vals(id);

        d1 = d * 3 * sqrt(n) * [6,5,4,2,3,7] ;
        d2 = d * 2 * sqrt(n) * [7,6,5,2,3,4,8,9,10] ;

        for ia = 1:na
            ang = ang_vals(ia);

            seed_dataset = BASE_SEED + 100000*rep_idx + 1000*id + 10*ia;
            rng(seed_dataset, 'twister');

            %% data generation
            v0 = random_orthonormal(n, r);
            v1 = random_orthonormal_orthogonal(r1 - r, v0);
            v_tmp = random_orthonormal_orthogonal(r2 - r, [v0 v1]);

            v2 = [ ...
                v_tmp(:, 1:(r1-r)) * sin(ang) + v1 * cos(ang), ...
                v_tmp(:, (r1-r+1):end) ...
            ];

            rot = random_orthonormal(r, r);

            X1 = u1 * diag(d1) * [v0, v1]'     + s1 * randn(p1, n);
            X2 = u2 * diag(d2) * [v0*rot, v2]' + s2 * randn(p2, n);
            
            true_j1 = u1 * diag([d1(1:r), zeros(1, r1-r)]) * [v0, v1]';
            true_j2 = u2 * diag([d2(1:r), zeros(1, r2-r)]) * [v0*rot, v2]';

            true_i1 = u1 * diag([zeros(1, r), d1(r+1:end)]) * [v0, v1]';
            true_i2 = u2 * diag([zeros(1, r), d2(r+1:end)]) * [v0*rot, v2]';

            true_s1 = true_j1 + true_i1;
            true_s2 = true_j2 + true_i2;

            relerr = @(A,B) norm(A-B,'fro')^2 / max(norm(A,'fro')^2, eps);

            %% methods
            for mm = 1:nm
                method = methods{mm};

                method_file = fullfile(raw_full_dir, ...
                    sprintf('rep_%03d_%s_result.mat', rep_idx, method));

                if exist(method_file, 'file') == 2
                    fprintf('[skip] rep %03d method %s already exists\n', rep_idx, method);
                    continue;
                end

                t0 = tic;

                rJ = NaN;
                rI1 = NaN;
                rI2 = NaN;
                r0_1 = NaN;
                r0_2 = NaN;

                J1_hat = zeros(size(X1));
                J2_hat = zeros(size(X2));
                I1_hat = zeros(size(X1));
                I2_hat = zeros(size(X2));
                S1_hat = zeros(size(X1));
                S2_hat = zeros(size(X2));

                err_msg = '';

                try
                    switch method

                        %% ==================  MSSAT ==================
                        case 'mssat'
                            [s1_hat2, r1_hat_j] = BEMA_combined(X1);
                            [s2_hat2, r2_hat_j] = BEMA_combined(X2);

                            r0_1 = r1_hat_j;
                            r0_2 = r2_hat_j;

                            X_list = {X1, X2};

                            out_mssat = mssat( ...
                                X_list, ...
                                'center', false, ...
                                'alpha', normcdf(-6), ...
                                'rank_method', 'given', ...
                                'sigma_method', 'given', ...
                                'r_init', [r1_hat_j, r2_hat_j], ...
                                'sigma_init', [s1_hat2, s2_hat2], ...
                                'test_mode', 'alpha', ...
                                'mode', 'test', ...
                                'alt_angle_deg', 5, ...
                                'verbose', false);

                            rJ = double(out_mssat.joint_rank);

                            if rJ > 0
                                Pj = out_mssat.P_joint;
                            else
                                Pj = zeros(n);
                            end

                            J1_hat = X1 * Pj;
                            J2_hat = X2 * Pj;

                            if isfield(out_mssat, 'blocks') && ~isempty(out_mssat.blocks)
                                if isfield(out_mssat.blocks{1}, 'P_indiv')
                                    Pind1 = out_mssat.blocks{1}.P_indiv;
                                    Pind2 = out_mssat.blocks{2}.P_indiv;
                                else
                                    Pind1 = out_mssat.blocks{1}.P_VAk - out_mssat.blocks{1}.P_VAk_given_VJ;
                                    Pind2 = out_mssat.blocks{2}.P_VAk - out_mssat.blocks{2}.P_VAk_given_VJ;
                                end

                                I1_hat = X1 * Pind1;
                                I2_hat = X2 * Pind2;
                            end

                            S1_hat = J1_hat + I1_hat;
                            S2_hat = J2_hat + I2_hat;

                            rI1 = max(r1_hat_j - rJ, 0);
                            rI2 = max(r2_hat_j - rJ, 0);

                        %% ================== AJIVE ==================
                        case 'ajive'
                            datablock = {X1, X2};
                            dataname = {'X','Y'};

                            [~, r1_hat] = BEMA_combined(X1);
                            [~, r2_hat] = BEMA_combined(X2);

                            r0_1 = r1_hat;
                            r0_2 = r2_hat;

                            vecr = [r1_hat, r2_hat];

                            paramstruct0 = struct( ...
                                'dataname', {dataname}, ...
                                'iplot', [0 0], ...
                                'ioutput', [0 0 0 0 0 0 1 1 0]);

                            outstruct0 = AJIVEMainMJ(datablock, vecr, paramstruct0);

                            rJ  = outstruct0.rjoint;
                            rI1 = rank(outstruct0.MatrixIndiv{1});
                            rI2 = rank(outstruct0.MatrixIndiv{2});

                            J1_hat = outstruct0.MatrixJoint{1};
                            J2_hat = outstruct0.MatrixJoint{2};

                            I1_hat = outstruct0.MatrixIndiv{1};
                            I2_hat = outstruct0.MatrixIndiv{2};

                            S1_hat = J1_hat + I1_hat;
                            S2_hat = J2_hat + I2_hat;

                        %% ================== DIVAS ==================
                        case 'divas'
                            data_sim = {X1, X2};

                            DIVASout = DJIVEMainJP(data_sim);

                            rJ = DIVASout.rjoint{1}(1);

                            keyVals = DIVASout.keyIdxMap.values;
                            keyVals = keyVals(:);

                            block_est = struct('inv',[],'joint',[]);

                            for ii = 1:2
                                vals = DIVASout.matBlocks{ii}.values;

                                if numel(vals) == 2
                                    block_est(ii).inv = vals{1};
                                    block_est(ii).joint = vals{2};
                                else
                                    is_single_ii = false;

                                    for kk = 1:numel(keyVals)
                                        try
                                            idxs = keyVals{kk}(:).';
                                            if numel(idxs) == 1 && idxs == ii
                                                is_single_ii = true;
                                                break;
                                            end
                                        catch
                                        end
                                    end

                                    if is_single_ii
                                        block_est(ii).inv = vals{1};
                                        block_est(ii).joint = zeros(size(vals{1}));
                                    else
                                        block_est(ii).inv = zeros(size(vals{1}));
                                        block_est(ii).joint = vals{1};
                                    end
                                end
                            end

                            rI1 = rank(block_est(1).inv);
                            rI2 = rank(block_est(2).inv);

                            J1_hat = block_est(1).joint;
                            J2_hat = block_est(2).joint;

                            I1_hat = block_est(1).inv;
                            I2_hat = block_est(2).inv;

                            S1_hat = J1_hat + I1_hat;
                            S2_hat = J2_hat + I2_hat;

                        otherwise
                            error('Unknown method: %s', method);
                    end

                catch ME
                    err_msg = ME.message;
                end

                dt = toc(t0);

                e_all_1   = relerr(true_s1, S1_hat);
                e_all_2   = relerr(true_s2, S2_hat);
                e_joint_1 = relerr(true_j1, J1_hat);
                e_joint_2 = relerr(true_j2, J2_hat);
                e_ind_1   = relerr(true_i1, I1_hat);
                e_ind_2   = relerr(true_i2, I2_hat);

                rep_result.(method).rJ(id, ia) = rJ;
                rep_result.(method).rI1(id, ia) = rI1;
                rep_result.(method).rI2(id, ia) = rI2;

                rep_result.(method).init_rank1(id, ia) = r0_1;
                rep_result.(method).init_rank2(id, ia) = r0_2;

                rep_result.(method).time(id, ia) = dt;

                rep_result.(method).err_overall(id, ia, :) = [e_all_1, e_all_2];
                rep_result.(method).err_joint(id, ia, :) = [e_joint_1, e_joint_2];
                rep_result.(method).err_ind(id, ia, :) = [e_ind_1, e_ind_2];

                rep_result.(method).error_msg{id, ia} = err_msg;
            end
        end
    end

    rep_result.elapsed = toc(rep_tic);

    % Save one file per method. This prevents later single-method runs
    % from overwriting previously computed results for other methods.
    for mm_save = 1:nm
        method_save = methods{mm_save};

        out_file = fullfile(raw_full_dir, ...
            sprintf('rep_%03d_%s_result.mat', rep_idx, method_save));

        if exist(out_file, 'file') == 2
            continue;
        end

        rep_method_result = struct();
        rep_method_result.rep_idx = rep_idx;
        rep_method_result.method = method_save;
        rep_method_result.d_vals = d_vals;
        rep_method_result.ang_vals = ang_vals;
        rep_method_result.result = rep_result.(method_save);
        rep_method_result.elapsed = rep_result.elapsed;

        parsave_rep_method_result(out_file, rep_method_result);
    end

    fprintf('[worker] rep %d / %d done, elapsed %.2fs\n', ...
        rep_idx, rep, rep_result.elapsed);
end

%% =========================================================
%% Aggregate results
%% =========================================================

fprintf('\nAggregating rep files...\n');

acc = struct();

for mm = 1:nm
    method = methods{mm};

    acc.(method).rJ = nan(rep, nd, na);
    acc.(method).rI1 = nan(rep, nd, na);
    acc.(method).rI2 = nan(rep, nd, na);

    acc.(method).time = nan(rep, nd, na);

    acc.(method).err_overall = nan(rep, nd, na, 2);
    acc.(method).err_joint = nan(rep, nd, na, 2);
    acc.(method).err_ind = nan(rep, nd, na, 2);

    acc.(method).init_rank1 = nan(rep, nd, na);
    acc.(method).init_rank2 = nan(rep, nd, na);
end

for rep_idx = 1:rep
    for mm = 1:nm
        method = methods{mm};

        in_file = fullfile(raw_full_dir, ...
            sprintf('rep_%03d_%s_result.mat', rep_idx, method));

        if exist(in_file, 'file') ~= 2
            warning('Missing file: %s', in_file);
            continue;
        end

        S = load(in_file, 'rep_method_result');
        rr = S.rep_method_result.result;

        acc.(method).rJ(rep_idx,:,:) = rr.rJ;
        acc.(method).rI1(rep_idx,:,:) = rr.rI1;
        acc.(method).rI2(rep_idx,:,:) = rr.rI2;

        acc.(method).time(rep_idx,:,:) = rr.time;

        acc.(method).err_overall(rep_idx,:,:,:) = rr.err_overall;
        acc.(method).err_joint(rep_idx,:,:,:) = rr.err_joint;
        acc.(method).err_ind(rep_idx,:,:,:) = rr.err_ind;

        acc.(method).init_rank1(rep_idx,:,:) = rr.init_rank1;
        acc.(method).init_rank2(rep_idx,:,:) = rr.init_rank2;
    end
end

%% =========================================================
%% Print final joint tables
%% =========================================================

for mm = 1:nm
    method = methods{mm};

    rJ_arr = acc.(method).rJ;
    time_arr = acc.(method).time;

    joint_mean = squeeze(mean(rJ_arr, 1, 'omitnan'));
    joint_prop = squeeze(mean(rJ_arr == true_r, 1, 'omitnan'));
    time_mean = squeeze(mean(time_arr, 1, 'omitnan'));

    C = cell(nd, 2*na);

    for id = 1:nd
        for ia = 1:na
            C{id, ia} = sprintf('%.2f(%.2f)', joint_mean(id, ia), joint_prop(id, ia));
            C{id, na+ia} = sprintf('%.4f', time_mean(id, ia));
        end
    end

    varNames = cell(1, 2*na);
    for ia = 1:na
        varNames{ia} = sprintf('J_ang%d', ia);
        varNames{na+ia} = sprintf('T_ang%d', ia);
    end

    rowNames = cellstr(num2str(d_vals(:), 'd=%.2f'));

    T = cell2table(C, 'VariableNames', varNames, 'RowNames', rowNames);

    fprintf('\n-------------------- %s FINAL JOINT --------------------\n', upper(method));
    disp([repmat({'joint'},1,na), repmat({'time'},1,na)]);
    disp(T);
end

save(fullfile(raw_full_dir, 'summary_parallel.mat'), ...
    'acc', 'methods', 'd_vals', 'ang_vals', 'rep', ...
    'true_r', 'r1_ind_true', 'r2_ind_true', '-v7.3');

fprintf('\nDONE. Summary saved:\n%s\n', fullfile(raw_full_dir, 'summary_parallel.mat'));

%% =========================================================
%% Print error tables separately
%% =========================================================

err_types = {'err_overall', 'err_joint', 'err_ind'};

for mm = 1:nm

    method = methods{mm};

    for ee = 1:length(err_types)

        err_name = err_types{ee};

        err_arr = acc.(method).(err_name);

        % size:
        % (rep, nd, na, 2)

        nMetric = size(err_arr, 4);

        for kk = 1:nMetric

            err_mean = squeeze(mean(err_arr(:,:,:,kk), 1, 'omitnan'));
            err_sd   = squeeze(std(err_arr(:,:,:,kk), 0, 1, 'omitnan'));

            C = cell(nd, na);

            for id = 1:nd
                for ia = 1:na

                    C{id, ia} = sprintf('%.6f(%.6f)', ...
                        err_mean(id, ia), ...
                        err_sd(id, ia));

                end
            end

            varNames = cell(1, na);

            for ia = 1:na
                varNames{ia} = sprintf('ang%d', ia);
            end

            rowNames = cellstr(num2str(d_vals(:), 'd=%.2f'));

            T = cell2table(C, ...
                'VariableNames', varNames, ...
                'RowNames', rowNames);

            fprintf('\n-------------------- %s %s metric%d --------------------\n', ...
                upper(method), upper(err_name), kk);

            disp(T);

        end
    end
end

%% =========================================================
%% local helper: parfor save
%% =========================================================

function parsave_rep_method_result(filename, rep_method_result)
    save(filename, 'rep_method_result', '-v7');
end