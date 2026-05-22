% Korean mortality data

var_keep = {'result_real', 'finalTable'};
vars = setdiff(who, var_keep);
clear(vars{:});
clc

addpath('~/MATLAB/MSSAT');

%% 
var_keep = {'time', 'realdata', 'jrank', 'result_real', 'finalTable'};
vars = setdiff(who, var_keep);
clear(vars{:});
clc

realdata = load('~/MATLAB/mydata/data_Male_Female_log.mat');

fields = fieldnames(realdata);  
numFields = numel(fields);      

dataCell = cell(1, numFields);
for i = 1:numFields
    tempData = realdata.(fields{i});
    
    if iscell(tempData) 
        dataCell{i} = cellfun(@str2double, tempData); 
    else
        dataCell{i} = tempData; 
    end
end

realdata = dataCell;

%% MSSAT

X_list = realdata;   % {X1, X2, X3}
K = numel(X_list);
fprintf('Running mssat on %d data blocks (Kor mortality)...\n', K);

tic
out_mssat = mssat(X_list, ...
    'center', true, ...
    'alpha', normcdf(-6), ...
    'rank_method', 'bema', ...
    'sigma_method', 'bema', ...
    'mode', 'test', ...
    'verbose', true);
time.mssat = toc;

if isfield(out_mssat, 'blocks') && ~isempty(out_mssat.blocks)
    r_individual = zeros(1, K);
    for k = 1:K
        try
            Vk_block = out_mssat.blocks{k}.P_VAk;
            r_individual(k) = rank(Vk_block);
        catch
            r_individual(k) = 0;
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


rJ = double(out_mssat.joint_rank);    % int32 → double 
r_individual = double(r_individual);
r_specific = max(r_individual - rJ, 0);    % block-specific ranks
r_specific = double(r_specific);

% ======================= results =======================
jrank.mssat = [rJ, r_specific];

disp('-------------------------------------------------------');
fprintf('>>> mssat finished.\n');
fprintf('Estimated joint rank: %d\n', rJ);
fprintf('Individual ranks (after joint removal): [%s]\n', num2str(r_specific));
disp('-------------------------------------------------------');

%% ===================== Paper-ready 2x2 (age-gradient, 20–39 color-emphasized) =====================
years = 1970:2023;      % 54
ages  = 0:100;          % 101

% ---- emphasize age band by COLOR remapping (not by fading others) ----
emphBand = [20 39];         
baseMap  = parula(256);      
tEmph    = [0.25 0.65];   

[cmapAgeEmph, t_of_age] = make_emphasis_colormap(ages, emphBand, baseMap, tEmph);

as_age_year = @(X) enforce_age_year(X, numel(ages), numel(years));

% ---- data: age x year ----
Xj_m = as_age_year(out_mssat.X_joint_hat{1});
Xj_f = as_age_year(out_mssat.X_joint_hat{2});
Xi_m = as_age_year(out_mssat.X_indiv_hat{1});
Xi_f = as_age_year(out_mssat.X_indiv_hat{2});

% ---- one figure, 2x2 ----
fig = figure('Color','w','Position',[60 60 1200 780]);
t = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile(t,1); plot_age_curves_colormap(ax1, years, ages, Xj_m, cmapAgeEmph, 'Male: Joint', true, emphBand);
ax2 = nexttile(t,2); plot_age_curves_colormap(ax2, years, ages, Xj_f, cmapAgeEmph, 'Female: Joint', true, emphBand);
ax3 = nexttile(t,3); plot_age_curves_colormap(ax3, years, ages, Xi_m, cmapAgeEmph, 'Male: Individual', false, emphBand);
ax4 = nexttile(t,4); plot_age_curves_colormap(ax4, years, ages, Xi_f, cmapAgeEmph, 'Female: Individual', false, emphBand);

% ---- shared colorbar (age) ----
colormap(fig, cmapAgeEmph);
cb = colorbar(ax4);
cb.Layout.Tile = 'east';
cb.Label.String = 'Age';
cb.Ticks = linspace(0,1,6);
cb.TickLabels = string(round(linspace(min(ages), max(ages), 6)));

title(t, 'Age-specific curves over years', ...
    'FontWeight','bold');


%% 

result_real.mortkor = combineRankAndTime(jrank, time);   % rank 1: joint rank, rank 2 and 3: individual rank

%% 

disp("result from mortality.kor");
disp(result_real.mortkor);


%% =========================================================
%% Korean mortality: original / joint / individual / residual
%% 4 x 2 heatmaps
%% =========================================================

K = numel(realdata);
if K ~= 2
    error('This code assumes exactly 2 blocks: Male/Female.');
end

block_names = ["Male", "Female"];
years = 1970:2023;
ages  = 0:100;

%% ---------- collect matrices ----------
X_original = cell(K,1);
X_joint    = cell(K,1);
X_indiv    = cell(K,1);
X_resid    = cell(K,1);

for k = 1:K
    X0 = double(realdata{k});

    % same centering convention as mssat(center=true)
    X_original{k} = bsxfun(@minus, X0, mean(X0, 2));

    X_joint{k} = double(out_mssat.X_joint_hat{k});
    X_indiv{k} = double(out_mssat.X_indiv_hat{k});

    X_resid{k} = X_original{k} - X_joint{k} - X_indiv{k};

    X_original{k}(~isfinite(X_original{k})) = 0;
    X_joint{k}(~isfinite(X_joint{k}))       = 0;
    X_indiv{k}(~isfinite(X_indiv{k}))       = 0;
    X_resid{k}(~isfinite(X_resid{k}))       = 0;
end

%% ---------- raw-value symmetric scale ----------
allvals = [
    X_original{1}(:);
    X_original{2}(:);
    X_joint{1}(:);
    X_joint{2}(:);
    X_indiv{1}(:);
    X_indiv{2}(:);
    X_resid{1}(:);
    X_resid{2}(:)
];

cmax_raw = prctile(abs(allvals), 99);
if isnan(cmax_raw) || cmax_raw <= 0
    cmax_raw = max(abs(allvals));
end
if isnan(cmax_raw) || cmax_raw <= 0
    cmax_raw = 1;
end

%% ---------- blue-white-red colormap ----------
m = 256;
n1 = floor(m/2);
n2 = m - n1;

blue_to_white = [ ...
    linspace(0,1,n1)', ...
    linspace(0,1,n1)', ...
    ones(n1,1) ...
];

white_to_red = [ ...
    ones(n2,1), ...
    linspace(1,0,n2)', ...
    linspace(1,0,n2)' ...
];

cmap_bwr = [blue_to_white; white_to_red];

%% ---------- ticks ----------
x_tick_pos = round(linspace(1, numel(years), 6));
x_tick_lab = string(years(x_tick_pos));

y_tick_pos = round(linspace(1, numel(ages), 6));
y_tick_lab = string(ages(y_tick_pos));

%% ---------- one combined figure: 2 x 4 ----------
fig_all = figure('Color','w','Position',[80 80 1550 650]);

t = tiledlayout(fig_all, 2, 4, ...
    'TileSpacing','compact', ...
    'Padding','compact');

plot_mats = {
    X_original{1}, "Male: Original";
    X_joint{1},    "Male: Joint";
    X_indiv{1},    "Male: Individual";
    X_resid{1},    "Male: Residual";

    X_original{2}, "Female: Original";
    X_joint{2},    "Female: Joint";
    X_indiv{2},    "Female: Individual";
    X_resid{2},    "Female: Residual"
};

ax_all = gobjects(8,1);

for i = 1:8
    ax = nexttile(t, i);
    ax_all(i) = ax;

    Xplot = plot_mats{i,1}; 

    imagesc(Xplot);
    axis tight;

    colormap(ax, cmap_bwr);
    clim(ax, [-cmax_raw cmax_raw]);

    xlabel('Year');
    ylabel('Age');

    xticks(x_tick_pos);
    xticklabels(x_tick_lab);

    yticks(y_tick_pos);
    yticklabels(y_tick_lab);

    title(plot_mats{i,2}, ...
        'Interpreter','none', ...
        'FontWeight','bold');

    set(gca, 'FontSize', 10);
end

%% ---------- two row-wise colorbars on the right ----------
cb1 = colorbar(ax_all(4));
cb2 = colorbar(ax_all(8));

cb1.Position = [0.965 0.533 0.006 0.405];  % [left bottom width height]
cb2.Position = [0.965 0.053 0.006 0.405];

title(t, 'Korean mortality', ...
    'FontSize', 14, ...
    'FontWeight', 'bold');

%% ---------- save ----------
mortality_mssat_heatmap.fig_all = fig_all;

