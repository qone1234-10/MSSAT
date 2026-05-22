function out = mssat(X_list, varargin)
% multi-block joint-direction scan
% X_list: 1xK cell, X_list{k} is p_k x n
%
% Options (Name-Value):
%   'center'       (logical) default true, if true do **row-centering** per block
%   'alpha'        (double)  default 0.05, one-sided test level (used when test_mode='alpha')
%   'verbose'      (logical) default true
%   'tol_proj'     (double)  default 1e-10
%   'rank_method'  'bema'|'mp'|'given'   (default 'bema')
%   'sigma_method' 'bema'|'mp'|'given'   (default 'bema')
%   'r_init'       (1xK int) required if rank_method='given'
%   'sigma_init'   (1xK dbl) required if sigma_method='given'
%
% Testing options:
%   'test_mode'     'alpha'|'alpha_seq'|'altmean' (default 'alpha')
%                   - 'alpha'     : use fixed alpha for every sequential test
%                   - 'alpha_seq' : use alpha_s = alpha / 2^s at the s-th sequential test
%                   - 'altmean'   : use alternative limit mean with fixed angle
%   'alt_angle_deg' (double) default 10
%                   - angle (deg) used for all pairwise true angles in Theorem 4 (altmean mode)
%   'mode'          'test'|'all' (default 'test')
%                   Kept for backward compatibility.
%                   The function now always scans all s=1..smax for test history.
%                   The reported joint rank always follows the original
%                   sequential rule: count only accepts before the first reject.
% ------------------------- Parse options -------------------------
p = inputParser;
addParameter(p, 'center', true, @(x)islogical(x)||ismember(x,[0,1]));
addParameter(p, 'alpha', 0.05, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p, 'verbose', true, @(x)islogical(x)||ismember(x,[0,1]));
addParameter(p, 'tol_proj', 1e-10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p, 'rank_method', 'bema', @(s)ischar(s)||isstring(s));
addParameter(p, 'sigma_method','bema', @(s)ischar(s)||isstring(s));
addParameter(p, 'r_init', [], @(x)isnumeric(x));
addParameter(p, 'sigma_init', [], @(x)isnumeric(x));

% New options
addParameter(p, 'test_mode', 'alpha', @(s)ischar(s)||isstring(s));
addParameter(p, 'alt_angle_deg', 10, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<180);
addParameter(p, 'mode', 'test', @(s)ischar(s)||isstring(s));

parse(p, varargin{:});
opt = p.Results;

rank_method  = validatestring(lower(string(opt.rank_method)), {'bema','mp','given'});
sigma_method = validatestring(lower(string(opt.sigma_method)), {'bema','mp','given'});
test_mode    = validatestring(lower(string(opt.test_mode)),    {'alpha','alpha_seq','altmean'});
mode         = validatestring(lower(string(opt.mode)),         {'test','all'});

K = numel(X_list);
assert(K>=2, 'K >= 2 required.');
n = size(X_list{1}, 2);
for k=1:K
    assert(size(X_list{k},2)==n, 'All blocks must share the same n (columns).');
end

% ------------------------- (0) Row-Centering -------------------------
if opt.center
    for k=1:K
        Xk = X_list{k};
        X_list{k} = bsxfun(@minus, Xk, mean(Xk,2));  % row-centering
    end
end

% ----------------- (1) Initial sigma^2 and ranks (once) -----------------
    function [sig2_mp, r_mp] = mp_sigma_rank(X)
        nk = size(X,2);
        pk = size(X,1);
        S  = (X.'*X)/nk;
        ev = eig((S+S.')/2);
        ev = max(ev,0);
        y  = pk/nk;
        sig2_mp = median(ev);                 % MP noise variance proxy
        edge    = sig2_mp * (1 + sqrt(y))^2;  % MP upper edge
        r_mp    = sum(ev > edge);
    end

    function r_mp_only = mp_rank(X)
        [~, r_mp_only] = mp_sigma_rank(X);
    end

    function b = est_block_init(X, kidx)
        sigma2_hat = [];
        rhat = [];

        % sigma/rank joint estimation logic
        if sigma_method=="mp"   && rank_method=="mp"
            [sigma2_hat, rhat] = mp_sigma_rank(X);

        elseif sigma_method=="mp" && rank_method~="mp"
            [sigma2_hat, ~] = mp_sigma_rank(X);

        elseif sigma_method~="mp" && rank_method=="mp"
            rhat = mp_rank(X);

        elseif sigma_method=="bema" && rank_method=="bema"
            [sigma2_hat, rhat] = BEMA_combined(X.', 0.2);

        elseif sigma_method=="bema" && rank_method~="bema"
            [sigma2_hat, ~] = BEMA_combined(X.', 0.2);

        elseif sigma_method~="bema" && rank_method=="bema"
            [~, rhat] = BEMA_combined(X.', 0.2);
        end

        % finalize sigma
        switch sigma_method
            case 'bema'
                sigma2_hat = max(sigma2_hat, 0);
            case 'mp'
                sigma2_hat = max(sigma2_hat, 0);
            case 'given'
                assert(numel(opt.sigma_init)==K, 'sigma_init must have length K.');
                sigma2_hat = max(opt.sigma_init(kidx), 0);
        end

        % finalize rank
        switch rank_method
            case 'bema'
                rhat = max(double(rhat), 0);
            case 'mp'
                rhat = max(double(rhat), 0);
            case 'given'
                assert(numel(opt.r_init)==K, 'r_init must have length K.');
                rhat = max(double(opt.r_init(kidx)), 0);
        end

        b.sigma2 = sigma2_hat;
        b.r      = rhat;
        b.y      = size(X,1)/size(X,2);   % p_k / n
    end

blks = cell(1,K);
for k=1:K
    blks{k} = est_block_init(X_list{k}, k);
end
rhat_vec = cellfun(@(b)b.r, blks);

if opt.verbose
    fprintf('Initial ranks (%s): %s\n', rank_method, ...
        strjoin(arrayfun(@(i) sprintf('r%d=%.0f',i,rhat_vec(i)), 1:K, 'UniformOutput',false), ', '));
    fprintf('Initial sigma^2 (%s): %s\n', sigma_method, ...
        strjoin(arrayfun(@(i) sprintf('%.4f', blks{i}.sigma2), 1:K, 'UniformOutput',false), ', '));
end

% ----------------- (2) Right-singular bases (row-subspaces) --------------
VA_list = cell(1,K);
for k=1:K
    rk = rhat_vec(k);
    if rk>0
        [~,~,V] = svd(X_list{k}, 'econ');
        VA_list{k} = V(:,1:rk);
    else
        VA_list{k} = zeros(n,0);
    end
end

% --------------------- (3) Plug-in helper functions ----------------------
inv_bgn = @(s2, sigma2, y) ...
    (sigma2/2) .* ( (s2./sigma2) - 1 - y + sqrt(max(((s2./sigma2) - 1 - y).^2 - 4*y, 0)) );
ell_hat = @(s, sigma2, y) ...
    max( sqrt(inv_bgn(s.^2, sigma2, y)) ./ (sqrt(n)*sqrt(sigma2)), 0);
a_y     = @(ell, y) ((ell.^4 - y) ./ (ell.^2 .* (ell.^2 + 1)));
theta_y = @(ell, y) ((ell.^4 + 2*y.*ell.^2 + y) ./ (ell.^3 .* (ell.^2 + 1).^2));
psi_y   = @(ell, y) ((ell.^6 - 3*y.*ell.^2 - 2*y) ./ (ell.^3 .* (ell.^2 + 1).^2));
V_E_y   = @(ell, y) ( 2./(ell.^4 - y) ) .* ( ...
        2*y.*(y+1).*theta_y(ell,y).^2 ...
      - (y.*(y-1).*(5*y+1))./(ell.*(ell.^2+1).^2).*theta_y(ell,y) ...
      + ((ell.^4 + y).*(ell.^2 + y).^2) ./ (ell.^3.*(ell.^2+1).^2) .* psi_y(ell,y) ...
      + 2*y.^2.*(y-1).^2 ./ (ell.^2.*(ell.^2+1).^4) );
V_y     = @(ell, y) 4*theta_y(ell,y).^2 + V_E_y(ell,y);

% ------------------ (4) One-pass scan over s = 1..smax -------------------
M = cell2mat(VA_list);
smax = min(rhat_vec);

% Prepare path logs (always allocate; will trim later if needed)
T2_path     = nan(max(smax,1),1);
crit_path   = nan(max(smax,1),1);
mu0_path    = nan(max(smax,1),1);
muAlt_path  = nan(max(smax,1),1);
V0_path     = nan(max(smax,1),1);
Z_path      = nan(max(smax,1),1);
accept_path = false(max(smax,1),1);

if smax == 0
    if opt.verbose
        fprintf('No singular vectors to test (smax=0).\n');
    end
    out = pack_output(0, zeros(n,0), zeros(n,n), [], {}, [], opt);

    out.V_block_joint  = repmat({zeros(n,0)}, 1, K);
    out.V_block_indiv  = VA_list;
    out.U_block_joint  = cell(1,K);
    out.U_block_indiv  = cell(1,K);
    for k=1:K
        pk = size(X_list{k},1);
        out.U_block_joint{k} = zeros(pk,0);
        out.U_block_indiv{k} = zeros(pk, size(VA_list{k},2));
    end

    out.test_path = struct('T2',[],'crit',[],'mu0',[],'mu_alt',[],'V0',[],'Z',[],'accept',[]);
    out.T2_path = [];
    out.crit_path = [];
    out.accept_path = [];
    return;
end

try
    [U_M,~,~] = svds(M, smax);
catch
    [Utmp,~,~] = svd(M, 'econ');
    U_M = Utmp(:,1:smax);
end

test_history = cell(0);
VJ_list = {};
accepted_s = [];

% Sequential-rank bookkeeping.
% We still evaluate all s=1,...,smax, but the estimated joint rank is
% determined only by the consecutive accepted prefix before the first reject.

prefix_alive = true;
first_reject_s = NaN;
rank_accept_path = false(max(smax,1),1);

for s = 1:smax
    vJ_hat = U_M(:, s);

    % =====================================================
    % PDF Algorithm 2 style block-specific lifted vectors
    %   w_hat     = Vhat_k' vJ_hat
    %   w_circ    = D_k^{-1} w_hat / ||D_k^{-1} w_hat||
    %   vtilde_k  = Vhat_k w_circ
    % =====================================================

    vAkJ    = cell(1,K);      % lifted vectors vtilde*_k
    m_hat   = zeros(K,1);     % m_hat(k)
    V_hat   = zeros(K,1);     % V_hat(k)

    a_tilde = zeros(K,1);     % stored as m_hat(k)^2 for compatibility
    V_tilde = zeros(K,1);     % stored as V_hat(k)

    tiny = 1e-12;

    for k = 1:K
        Vk = VA_list{k};
        rk = size(Vk, 2);

        if isempty(Vk) || rk == 0
            vAkJ{k} = zeros(n,1);
            m_hat(k) = 0;
            V_hat(k) = 0;
            a_tilde(k) = 0;
            V_tilde(k) = 0;
            continue;
        end

        Xk   = X_list{k};
        sig2 = blks{k}.sigma2;
        yk   = blks{k}.y;

        % candidate direction coordinates in block-k row subspace
        w_hat = Vk' * vJ_hat/norm(Vk' * vJ_hat);       % r_k x 1

        % singular values along Vk directions
        sdir_i = vecnorm(Xk * Vk, 2, 1)';  % r_k x 1

        % estimated spike strength and asymptotic quantities
        ell_i = ell_hat(sdir_i, sig2, yk);

        a_i = max(a_y(ell_i, yk), 0);
        V_i = max(V_y(ell_i, yk), 0);

        % m_y(lambda) is represented by sqrt(a_i)
        m_i = sqrt(max(a_i, 0));

        % numerical guard
        valid = isfinite(m_i) & isfinite(V_i) & (m_i > tiny);

        if ~any(valid)
            vAkJ{k} = zeros(n,1);
            m_hat(k) = 0;
            V_hat(k) = 0;
            a_tilde(k) = 0;
            V_tilde(k) = 0;
            continue;
        end

        m_i_safe = m_i;
        m_i_safe(~valid) = Inf;

        % D^{-1} w_hat
        Dinv_w = w_hat ./ m_i_safe;
        Dinv_w(~isfinite(Dinv_w)) = 0;

        norm_Dinv_w = norm(Dinv_w);

        if norm_Dinv_w <= tiny
            vAkJ{k} = zeros(n,1);
            m_hat(k) = 0;
            V_hat(k) = 0;
            a_tilde(k) = 0;
            V_tilde(k) = 0;
            continue;
        end

        % w_circ = D^{-1}w / ||D^{-1}w||
        w_circ = Dinv_w / norm_Dinv_w;

        % lifted vector vtilde*_k
        vAkJ{k} = Vk * w_circ;

        % m_hat(k) = sum_l w_l w_circ_l m_y(lambda_l)
        % m_hat(k) = sum(w_hat .* w_circ .* m_i);
        m_hat(k) = 1/norm_Dinv_w;
        

        % V_hat(k) = sum_l w_l^2 w_circ_l^2 V_y(lambda_l)
        V_hat(k) = sum((w_hat.^2) .* (w_circ.^2) .* V_i);

        % compatibility with previous variable names
        a_tilde(k) = m_hat(k)^2;
        V_tilde(k) = V_hat(k);
    end

    % =====================================================
    % Test statistic T^2 = || sum_k vtilde*_k ||^2
    % =====================================================

    Vsum = zeros(n,1);
    for k = 1:K
        Vsum = Vsum + vAkJ{k};
    end
    T2 = sum(Vsum.^2);

    % =====================================================
    % Plug-in null mean and variance
    %
    % mu0 = K + 2 sum_{i<j} m_i m_j
    %
    % V0 = 4 sum_k V_k (sum_{l != k} m_l)^2
    %      + 4 sum_{i<j} (1 - m_i^2)(1 - m_j^2)
    % =====================================================

    mK = min(max(m_hat, 0), 1);

    pair_sum = 0;
    if K >= 2
        for i = 1:K-1
            for j = i+1:K
                pair_sum = pair_sum + mK(i) * mK(j);
            end
        end
        mu0 = K + 2 * pair_sum;
    else
        mu0 = K;
    end

    cosphi = cosd(opt.alt_angle_deg);
    mu_alt = K + 2 * cosphi * pair_sum;

    V0 = 0;

    for i = 1:K
        others = [1:i-1, i+1:K];
        s_others = sum(mK(others));
        V0 = V0 + 4 * V_hat(i) * (s_others^2);
    end

    for i = 1:K-1
        for j = i+1:K
            V0 = V0 + 4 * (1 - mK(i)^2) * (1 - mK(j)^2);
        end
    end

    V0 = max(V0, 0);

    % =====================================================
    % Testing rule
    %
    % Under H0, Z = sqrt(n/V0) * (T2 - mu0) ~ N(0,1).
    % Misalignment shifts Z left.
    % Reject alignment if Z < z_alpha.
    % Accept joint direction if Z >= z_alpha.
    % =====================================================

    Zn   = NaN;
    crit = NaN;
    pval = NaN;
    alpha_s_used = NaN;

    switch test_mode
        case 'alpha'
            if V0 > tiny
                Zn = sqrt(n / V0) * (T2 - mu0);

                alpha_s = opt.alpha;
                alpha_s_used = alpha_s;

                crit = norminv(alpha_s, 0, 1);
                pval = normcdf(Zn, 0, 1);
            end

            accept_joint = (~isnan(Zn)) && (Zn >= crit);

        case 'alpha_seq'
            if V0 > tiny
                Zn = sqrt(n / V0) * (T2 - mu0);

                alpha_s = opt.alpha / (2^s);
                alpha_s_used = alpha_s;

                crit = norminv(alpha_s, 0, 1);
                pval = normcdf(Zn, 0, 1);
            end

            accept_joint = (~isnan(Zn)) && (Zn >= crit);

        case 'altmean'
            crit = mu_alt;

            if V0 > tiny
                Zn = sqrt(n / V0) * (T2 - mu0);
            end

            pval = NaN;

            accept_joint = T2 > crit;
    end

    % =====================================================
    % Record path
    % =====================================================

    T2_path(s)     = T2;
    crit_path(s)   = crit;
    mu0_path(s)    = mu0;
    muAlt_path(s)  = mu_alt;
    V0_path(s)     = V0;
    Z_path(s)      = Zn;
    accept_path(s) = logical(accept_joint);

    entry = struct( ...
        's', s, ...
        'T2', T2, ...
        'mu', mu0, ...
        'mu_alt', mu_alt, ...
        'V', V0, ...
        'Z', Zn, ...
        'p_value', pval, ...
        'crit', crit, ...
        'alpha_s', alpha_s_used, ...
        'm_hat', m_hat, ...
        'a_tilde', a_tilde, ...
        'V_tilde', V_tilde, ...
        'accept', logical(accept_joint) ...
        );
    test_history{end+1} = entry; %#ok<AGROW>

    if opt.verbose
        fprintf('Scan s=%d: T2=%.3f, Z=%.3f, crit=%.3f, mode=%s/%s, accept=%s\n', ...
                s, T2, Zn, crit, test_mode, mode, string(accept_joint));
    end

    if opt.verbose
        fprintf('Scan s=%d: T2=%.3f, Z=%.3f, crit=%.3f, mode=%s/%s, accept=%s\n', ...
            s, T2, Zn, crit, test_mode, mode, string(accept_joint));
    end

% =====================================================
% Sequential-rank rule with full path logging
%
% We do NOT stop the loop at the first rejection, because we want
% test_history and test_path for all s = 1,...,smax.
%
% However, the joint rank is computed exactly as in the original
% sequential algorithm: only the accepted directions before the first
% rejection are counted.
% =====================================================

    if prefix_alive && accept_joint
        VJ_list{end+1} = vJ_hat; %#ok<AGROW>
        accepted_s = [accepted_s, s]; %#ok<AGROW>
        rank_accept_path(s) = true;

    elseif prefix_alive && ~accept_joint
        first_reject_s = s;
        prefix_alive = false;
    end
end

% Trim paths to evaluated range
last_s = numel(test_history);
T2_path     = T2_path(1:last_s);
crit_path   = crit_path(1:last_s);
mu0_path    = mu0_path(1:last_s);
muAlt_path  = muAlt_path(1:last_s);
V0_path     = V0_path(1:last_s);
Z_path      = Z_path(1:last_s);
accept_path = accept_path(1:last_s);
rank_accept_path = rank_accept_path(1:last_s);

% ----------------------- (5) Pack results + block-wise V/U ------------------------
rJ_hat = numel(VJ_list);

if rJ_hat>0
    VJ_raw = cell2mat(VJ_list);
    [Uj,~,~] = svd(VJ_raw, 'econ');
    V_joint = Uj(:,1:rJ_hat);
    P_joint  = V_joint * V_joint.';

    rec = cell(1,K);
    V_block_joint = cell(1,K);
    V_block_indiv = cell(1,K);
    U_block_joint = cell(1,K);
    U_block_indiv = cell(1,K);

    for k=1:K
        Vk0 = VA_list{k};
        pk  = size(X_list{k},1);

        if ~isempty(Vk0) && size(Vk0,2)>0
            PVAk = Vk0 * Vk0.';
        else
            PVAk = zeros(n,n);
        end

        % block k joint part (row-side)
        if any(PVAk(:)~=0)
            Ak = PVAk * V_joint;
            [UA,DA,~] = svd(Ak, 'econ');
            r_eff = sum(diag(DA) > opt.tol_proj);
            if r_eff>0
                Vj_block = UA(:,1:r_eff);
                Pcond    = Vj_block * Vj_block.';
            else
                Vj_block = zeros(n,0);
                Pcond    = zeros(n,n);
            end
        else
            Vj_block = zeros(n,0);
            Pcond    = zeros(n,n);
        end

        % individual projector (row-side)
        Pind = PVAk - Pcond;

        % individual row basis
        if any(Pind(:) ~= 0)
            [Ui,Si,~] = svd(Pind, 'econ');
            r_ind = sum(diag(Si) > opt.tol_proj);
            if r_ind > 0
                Vind = Ui(:,1:r_ind);
            else
                Vind = zeros(n,0);
            end
        else
            Vind = zeros(n,0);
        end

        V_block_joint{k} = Vj_block;
        V_block_indiv{k} = Vind;

        % joint U for block k: col space of X_k * V_joint_block
        if ~isempty(Vj_block) && size(Vj_block,2)>0
            Xj_k = X_list{k} * Vj_block;     % p_k x r_{Jk}
            [Ujb,~,~] = svd(Xj_k, 'econ');
            U_block_joint{k} = Ujb(:,1:size(Vj_block,2));
        else
            U_block_joint{k} = zeros(pk,0);
        end

        % individual U for block k: col space of X_k * V_indiv_block
        if ~isempty(Vind) && size(Vind,2)>0
            Xi_k = X_list{k} * Vind;
            [Uib,~,~] = svd(Xi_k, 'econ');
            U_block_indiv{k} = Uib(:,1:size(Vind,2));
        else
            U_block_indiv{k} = zeros(pk,0);
        end

        rec{k} = struct( ...
            'P_VAk',           PVAk, ...
            'P_VAk_given_VJ',  Pcond, ...
            'P_indiv',         Pind, ...
            'V_joint_block',   Vj_block, ...
            'V_indiv',         Vind, ...
            'U_joint_block',   U_block_joint{k}, ...
            'U_indiv',         U_block_indiv{k});
    end

    out = pack_output(rJ_hat, V_joint, P_joint, accepted_s, test_history, rec, opt);
    out.V_block_joint = V_block_joint;
    out.V_block_indiv = V_block_indiv;
    out.U_block_joint = U_block_joint;
    out.U_block_indiv = U_block_indiv;

else

    V_joint       = zeros(n,0);
    P_joint       = zeros(n,n);
    rec           = cell(1,K);
    V_block_joint = cell(1,K);
    V_block_indiv = cell(1,K);
    U_block_joint = cell(1,K);
    U_block_indiv = cell(1,K);

    for k=1:K
        Vk0 = VA_list{k};
        pk  = size(X_list{k},1);

        if ~isempty(Vk0) && size(Vk0,2)>0
            PVAk = Vk0 * Vk0.';
        else
            PVAk = zeros(n,n);
        end

        Pcond    = zeros(n,n);
        Pind     = PVAk;
        Vj_block = zeros(n,0);
        Vind     = Vk0;

        V_block_joint{k} = Vj_block;
        V_block_indiv{k} = Vind;

        U_block_joint{k} = zeros(pk,0);
        if ~isempty(Vk0) && size(Vk0,2)>0
            Xi_k = X_list{k} * Vk0;
            [Uib,~,~] = svd(Xi_k, 'econ');
            U_block_indiv{k} = Uib(:,1:size(Vk0,2));
        else
            U_block_indiv{k} = zeros(pk,0);
        end

        rec{k} = struct( ...
            'P_VAk',           PVAk, ...
            'P_VAk_given_VJ',  Pcond, ...
            'P_indiv',         Pind, ...
            'V_joint_block',   Vj_block, ...
            'V_indiv',         Vind, ...
            'U_joint_block',   U_block_joint{k}, ...
            'U_indiv',         U_block_indiv{k});
    end

    out = pack_output(0, V_joint, P_joint, accepted_s, test_history, rec, opt);
    out.V_block_joint = V_block_joint;
    out.V_block_indiv = V_block_indiv;
    out.U_block_joint = U_block_joint;
    out.U_block_indiv = U_block_indiv;
end

% ===================== (6) Reconstruct joint/individual matrices =====================
% Reconstruct in data space (p_k x n) using right-side projectors in score space.

X_joint_hat  = cell(1,K);
X_indiv_hat  = cell(1,K);
X_signal_hat = cell(1,K);

for k = 1:K
    Xk = X_list{k};

    % joint reconstruction: Xk * P_joint (works even if rJ=0 because P_joint could be zeros)
    if isfield(out,'P_joint') && ~isempty(out.P_joint)
        X_joint_hat{k} = Xk * out.P_joint;
    else
        X_joint_hat{k} = zeros(size(Xk));
    end

    % individual reconstruction: Xk * P_indiv (block-specific)
    if isfield(out,'blocks') && ~isempty(out.blocks) && numel(out.blocks)>=k && isfield(out.blocks{k},'P_indiv')
        Pind = out.blocks{k}.P_indiv;
        if ~isempty(Pind)
            X_indiv_hat{k} = Xk * Pind;
        else
            X_indiv_hat{k} = zeros(size(Xk));
        end
    else
        % fallback: if blocks not available, use (I - P_joint)
        n = size(Xk,2);
        if isfield(out,'P_joint') && ~isempty(out.P_joint)
            X_indiv_hat{k} = Xk * (eye(n) - out.P_joint);
        else
            X_indiv_hat{k} = zeros(size(Xk));
        end
    end

    % estimated signal part (joint + individual)
    X_signal_hat{k} = X_joint_hat{k} + X_indiv_hat{k};
end

out.X_joint_hat  = X_joint_hat;
out.X_indiv_hat  = X_indiv_hat;
out.X_signal_hat = X_signal_hat;

% Attach path vectors for easy access / plotting
out.test_path = struct( ...
    'T2',     T2_path(:), ...
    'crit',   crit_path(:), ...
    'mu0',    mu0_path(:), ...
    'mu_alt', muAlt_path(:), ...
    'V0',     V0_path(:), ...
    'Z',      Z_path(:), ...
    'accept', accept_path(:) , ...
    'rank_accept', rank_accept_path(:) );

out.T2_path     = T2_path(:);
out.crit_path   = crit_path(:);
out.accept_path = accept_path(:);
out.first_reject_s   = first_reject_s;
out.rank_accept_path = rank_accept_path(:);

end

% ======================= Helper: pack_output =======================
function out = pack_output(rJ, VJ, PJ, acc_s, hist, rec, opt)
out = struct();
out.joint_rank   = double(rJ);
out.V_joint      = VJ;
out.P_joint      = PJ;
out.accepted_s   = double(acc_s(:)).';
out.test_history = hist;
out.blocks       = rec;
out.options      = struct('center',logical(opt.center), ...
                          'alpha',opt.alpha, ...
                          'tol_proj',opt.tol_proj, ...
                          'rank_method',char(opt.rank_method), ...
                          'sigma_method',char(opt.sigma_method), ...
                          'r_init',opt.r_init, ...
                          'sigma_init',opt.sigma_init, ...
                          'test_mode',char(opt.test_mode), ...
                          'alt_angle_deg',opt.alt_angle_deg, ...
                          'mode',char(opt.mode));
end