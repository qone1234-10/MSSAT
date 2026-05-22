%% ================== Setup ==================
clear; clc;

addpath('~/MATLAB/MSSAT');
addpath('~/MATLAB/MyJIVE');
addpath('~/MATLAB/AJIVE_Project/AJIVECode');
addpath('~/MATLAB/Data-Integration-Via-Analysis-of-Subspaces/DJIVECode');
addpath('~/MATLAB/cvx');

%% ================== USER OPTIONS ==================
% 'ajive', 'divas', 'mssat'

% methods = {'ajive', 'divas', 'mssat'};
methods = {'mssat'};

silence_method_output = true;     
print_per_method_progress = true; 
print_every_rep_table = true;     
print_final_tables = true;        

%% ================== Common parameters (fixed d) ==================
d  = 1;              % signal strength scale
r1 = 6; r2 = 9; r = 3;

s1 = 1; s2 = 1;

true_r      = r;
r1_ind_true = r1 - true_r;
r2_ind_true = r2 - true_r;

%% ================== Grid: n and aspect ratio ==================
n_vals  = [100 200 400];
ar_vals = [0.25 0.5 1 2 4];
rep     = 100;

nn  = numel(n_vals);
nar = numel(ar_vals);

fprintf('\n================= Simulation Grid =================\n');
fprintf('Methods order = %s\n', strjoin(methods, ', '));
fprintf('rep     = %d\n', rep);
fprintf('n_vals  = [%s]\n', num2str(n_vals));
fprintf('ar_vals = [%s]\n', num2str(ar_vals));
fprintf('Total datasets = nar*nn*rep = %d*%d*%d = %d\n', nar, nn, rep, nar*nn*rep);
fprintf('===================================================\n\n');

%% ================== Accumulators across reps ==================
% acc.(method).joint_sum(ar_idx, n_idx), joint_correct(...)
% acc.(method).ind_sum(ar_idx, n_idx, block), ind_correct(...)
% acc.(method).time_sum(ar_idx, n_idx)
acc = struct();
for m = 1:numel(methods)
    method = methods{m};

    acc.(method).joint_sum     = zeros(nar, nn);
    acc.(method).joint_correct = zeros(nar, nn);

    acc.(method).ind_sum       = zeros(nar, nn, 2);
    acc.(method).ind_correct   = zeros(nar, nn, 2);

    acc.(method).time_sum      = zeros(nar, nn);
end

overall_tic = tic;

%% ================== MAIN LOOP: rep OUTERMOST ==================
for rep_idx = 1:rep
    rep_tic = tic;

    fprintf('\n===================================================\n');
    fprintf('>>> REP %d / %d START (elapsed total %.1fs)\n', rep_idx, rep, toc(overall_tic));
    fprintf('===================================================\n');

    rng(0);

    for ia = 1:nar
        ar = ar_vals(ia);

        for in = 1:nn
            n = n_vals(in);

            % p1, p2 결정
            p1 = max(round(ar * n), r1 + 1);
            p2 = max(round(p1 * 4/5), r2 + 1);

            % singular values
            d1 = d * 3/2 * sqrt(n) * [6,5,4,2,3,7] ;
            d2 = d * sqrt(n) * [7,6,5,2,3,4,8,9,10] ;
            
            % subspaces
            v0    = random_orthonormal(n, r);
            v_tmp = random_orthonormal_orthogonal(r1 + r2 - 2*r, v0);
            v1 = v_tmp(:, 1:(r1 - r));
            v2 = v_tmp(:, (r1 - r + 1):end);

            % left singular vectors depend on p1/p2
            u1 = random_orthonormal(p1, r1);
            u2 = random_orthonormal(p2, r2);

            % data
            rot = random_orthonormal(r, r);
            X1  = u1 * diag(d1) * [v0, v1]'     + s1 * randn(p1, n);
            X2  = u2 * diag(d2) * [v0*rot, v2]' + s2 * randn(p2, n);

            % ======== METHODS (fixed order) ========
            for mm = 1:numel(methods)
                method = methods{mm};
                t0 = tic;

                if print_per_method_progress
                    fprintf('rep %3d/%3d | method=%-7s | ar=%.2f (%d/%d) | n=%3d (%d/%d) ... ', ...
                        rep_idx, rep, method, ar, ia, nar, n, in, nn);
                end

                switch method
                    case 'mssat'  
                        [s1_hat2, r1_hat_j] = BEMA_combined(X1);
                        [s2_hat2, r2_hat_j] = BEMA_combined(X2);

                        X_list = {X1, X2};

                        if silence_method_output
                            [~, out_mssat] = evalc(['mssat(X_list, ' ...
                                '''center'', false, ''alpha'', normcdf(-6), ' ...
                                '''rank_method'', ''given'', ''sigma_method'', ''given'', ' ...
                                '''r_init'', [r1_hat_j, r2_hat_j], ''sigma_init'', [s1_hat2, s2_hat2], ' ...
                                '''test_mode'', ''alpha'', ''alt_angle_deg'', 5, ''verbose'', false);']);
                        else
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
                        end

                        rJ = double(out_mssat.joint_rank);

                        if rJ > 0 && isfield(out_mssat,'blocks') && ~isempty(out_mssat.blocks)
                            Pind1 = out_mssat.blocks{1}.P_indiv;
                            Pind2 = out_mssat.blocks{2}.P_indiv;
                            rI1 = rank(Pind1);
                            rI2 = rank(Pind2);
                        else
                            rI1 = r1_hat_j;
                            rI2 = r2_hat_j;
                        end

                    case 'ajive'
                        datablock = {X1, X2};
                        dataname  = {'X','Y'};

                        [~, r1_hat] = BEMA_combined(X1);
                        [~, r2_hat] = BEMA_combined(X2);
                        vecr = [r1_hat, r2_hat];

                        paramstruct0 = struct('dataname', {dataname}, ...
                                              'iplot', [0 0], ...
                                              'ioutput', [0 0 0 0 0 0 1 1 0]);

                        if silence_method_output
                            [~, outstruct0] = evalc('AJIVEMainMJ(datablock, vecr, paramstruct0);');
                        else
                            outstruct0 = AJIVEMainMJ(datablock, vecr, paramstruct0);
                        end

                        rJ  = outstruct0.rjoint;
                        rI1 = rank(outstruct0.MatrixIndiv{1});
                        rI2 = rank(outstruct0.MatrixIndiv{2});

                    case 'divas'
                        data_sim = {X1, X2};

                        if silence_method_output
                            [~, DIVASout] = evalc('DJIVEMainJP(data_sim);');
                        else
                            DIVASout = DJIVEMainJP(data_sim);
                        end

                        rJ = DIVASout.rjoint{1}(1);

                        % robust indiv parsing (single-block collection detection)
                        keyVals = DIVASout.keyIdxMap.values; keyVals = keyVals(:);

                        block_est = struct('inv',[],'joint',[]);
                        for ii = 1:2
                            vals = DIVASout.matBlocks{ii}.values;

                            if numel(vals) == 2
                                block_est(ii).inv   = vals{1};
                                block_est(ii).joint = vals{2};
                            else
                                is_single_ii = false;
                                for kk = 1:numel(keyVals)
                                    try
                                        idxs = keyVals{kk}(:).';
                                        if numel(idxs)==1 && idxs==ii
                                            is_single_ii = true; break;
                                        end
                                    catch
                                    end
                                end
                                if is_single_ii
                                    block_est(ii).inv   = vals{1};
                                    block_est(ii).joint = zeros(size(vals{1}));
                                else
                                    block_est(ii).inv   = zeros(size(vals{1}));
                                    block_est(ii).joint = vals{1};
                                end
                            end
                        end
                        rI1 = rank(block_est(1).inv);
                        rI2 = rank(block_est(2).inv);

                    otherwise
                        error('Unknown method: %s', method);
                end

                dt = toc(t0);

                % accumulate across reps
                acc.(method).joint_sum(ia,in) = acc.(method).joint_sum(ia,in) + rJ;
                acc.(method).joint_correct(ia,in) = acc.(method).joint_correct(ia,in) + double(rJ == true_r);

                acc.(method).ind_sum(ia,in,1) = acc.(method).ind_sum(ia,in,1) + rI1;
                acc.(method).ind_sum(ia,in,2) = acc.(method).ind_sum(ia,in,2) + rI2;
                acc.(method).ind_correct(ia,in,1) = acc.(method).ind_correct(ia,in,1) + double(rI1 == r1_ind_true);
                acc.(method).ind_correct(ia,in,2) = acc.(method).ind_correct(ia,in,2) + double(rI2 == r2_ind_true);

                acc.(method).time_sum(ia,in) = acc.(method).time_sum(ia,in) + dt;

                if print_per_method_progress
                    fprintf('dt=%.3fs | rJ=%d | rep_elapsed=%.1fs\n', dt, rJ, toc(rep_tic));
                end
            end % methods

        end % n loop
    end % ar loop

    fprintf('\n>>> REP %d / %d END (elapsed this rep %.2fs, elapsed total %.1fs)\n', ...
        rep_idx, rep, toc(rep_tic), toc(overall_tic));

    %% ================== REP-END: JOINT-only interim tables ==================
    if print_every_rep_table
        fprintf('\n==================== Interim JOINT tables (reps=1..%d) ====================\n', rep_idx);

        for mm = 1:numel(methods)
            method = methods{mm};

            joint_mean = acc.(method).joint_sum     / rep_idx;   % nar x nn
            joint_prop = acc.(method).joint_correct / rep_idx;   % nar x nn
            time_mean  = acc.(method).time_sum      / rep_idx;   % nar x nn

            C = cell(nar, 2*nn);
            for ia2 = 1:nar
                for in2 = 1:nn
                    % 요청: 평균 rank(정확비율)
                    C{ia2, in2}      = sprintf('%.2f(%.2f)', joint_mean(ia2,in2), joint_prop(ia2,in2));
                    C{ia2, nn+in2}   = sprintf('%.4f', time_mean(ia2,in2));
                end
            end

            varNames = cell(1, 2*nn);
            for in2 = 1:nn
                varNames{in2}      = sprintf('J_n%d', n_vals(in2));
                varNames{nn+in2}   = sprintf('T_n%d', n_vals(in2));
            end
            rowNames = cellstr(num2str(ar_vals(:), 'ar=%.2f'));

            T = cell2table(C, 'VariableNames', varNames, 'RowNames', rowNames);

            fprintf('\n-------------------- %s (JOINT only, reps=1..%d) --------------------\n', upper(method), rep_idx);
            disp([repmat({'joint'},1,nn), repmat({'time'},1,nn)]);
            disp(T);
        end
    end
end % rep loop

%% ================== FINAL tables (include individuals) ==================
if print_final_tables
    fprintf('\n\n==================== FINAL tables (reps=1..%d) ====================\n', rep);

    for mm = 1:numel(methods)
        method = methods{mm};

        joint_mean = acc.(method).joint_sum     / rep;
        joint_prop = acc.(method).joint_correct / rep;

        ind_mean   = acc.(method).ind_sum       / rep;
        ind_prop   = acc.(method).ind_correct   / rep;

        time_mean  = acc.(method).time_sum      / rep;

        C = cell(nar, 4*nn);
        for ia = 1:nar
            for in = 1:nn
                C{ia, in}           = sprintf('%.2f(%.2f)', joint_mean(ia,in), joint_prop(ia,in));
                C{ia, nn+in}        = sprintf('%.2f(%.2f)', ind_mean(ia,in,1), ind_prop(ia,in,1));
                C{ia, 2*nn+in}      = sprintf('%.2f(%.2f)', ind_mean(ia,in,2), ind_prop(ia,in,2));
                C{ia, 3*nn+in}      = sprintf('%.4f', time_mean(ia,in));
            end
        end

        varNames = cell(1, 4*nn);
        for in = 1:nn
            varNames{in}          = sprintf('J_n%d',  n_vals(in));
            varNames{nn+in}       = sprintf('I1_n%d', n_vals(in));
            varNames{2*nn+in}     = sprintf('I2_n%d', n_vals(in));
            varNames{3*nn+in}     = sprintf('T_n%d',  n_vals(in));
        end
        rowNames = cellstr(num2str(ar_vals(:), 'ar=%.2f'));

        Tfinal = cell2table(C, 'VariableNames', varNames, 'RowNames', rowNames);

        fprintf('\n-------------------- %s (FINAL, reps=1..%d) --------------------\n', upper(method), rep);
        disp([repmat({'joint'},1,nn), repmat({'ind1'},1,nn), repmat({'ind2'},1,nn), repmat({'time'},1,nn)]);
        disp(Tfinal);
    end

    fprintf('\n================= ALL DONE =================\n');
    fprintf('Total elapsed time = %.2f sec\n', toc(overall_tic));
    fprintf('===========================================\n');
end
