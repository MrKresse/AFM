function results = batch_processing_mfjc(path_processed, force_dir, ...
        L_cont0, T0, l_kuhn0, k_seg0, v0, varargin)
% BATCH_PROCESSING_MFJC
% Loads rupture intervals from rupture_info.mat, fits only *_rupture.mat curves
% with fit_mFJC, and saves results in UI-compatible format (fitResults.mat).
%
% Inputs:
%   path_processed : folder with processed data and rupture_info.mat
%   force_dir      : struct array from dir(...), used for base names
%   L_cont0, T0, l_kuhn0, k_seg0, v0 : initial guesses / constants
%
% Optional name-value parameters:
%   'cache_file'        [fullfile(path_processed,'rupture_info.mat')]
%   'SavePlots'         [true]
%   'OutDir'            [fullfile(path_processed,'mFJC_fits')]
%   'ReplaceExisting'   [true]
%   'Quiet'             [false]
%
% Returns:
%   results : table summary of fitted rupture events

% ---------------- Options ----------------
p = inputParser;
p.addParameter('cache_file', fullfile(path_processed,'rupture_info.mat'));
p.addParameter('SavePlots', true, @(x)islogical(x));
p.addParameter('OutDir', fullfile(path_processed,'mFJC_fits'), @ischar);
p.addParameter('ReplaceExisting', true, @(x)islogical(x));
p.addParameter('Quiet', false, @(x)islogical(x));
p.parse(varargin{:});
opt = p.Results;

if opt.SavePlots && ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

% ---------------- Load rupture cache ----------------
cache_file = opt.cache_file;
if ~isfile(cache_file)
    error('batch_processing_mfjc:cache_missing', ...
        'Cache file not found: %s. Run batch_collect_ruptures first.', cache_file);
end
S = load(cache_file);
if ~isfield(S,'rupture_info') || isempty(S.rupture_info)
    error('batch_processing_mfjc:cache_empty', ...
        'Cache file %s does not contain rupture_info or it is empty.', cache_file);
end
rupture_info = S.rupture_info;

% build base_str -> [z0_left, z_pre]
info_map = containers.Map('KeyType','char','ValueType','any');
for k = 1:numel(rupture_info)
    if isfield(rupture_info(k),'base_str') && ~isempty(rupture_info(k).base_str)
        if isfield(rupture_info(k),'z0_left') && isfield(rupture_info(k),'z_pre')
            info_map(rupture_info(k).base_str) = [rupture_info(k).z0_left, rupture_info(k).z_pre];
        end
    end
end

% ---------------- Restrict to rupture files only ----------------
all_files = dir(fullfile(path_processed, '*_rupture.mat'));
if isempty(all_files)
    error('No *_rupture.mat files found in %s', path_processed);
end
rupture_names = {all_files.name};
% strip suffix to get the base names
bases_with_rupture = erase(rupture_names, '_retract_scaled.mat_rupture.mat');

% ---------------- Prepare outputs ----------------
fit_path = fullfile(path_processed, 'fitResults.mat');
fitResults = [];
if exist(fit_path,'file')
    Tfit = load(fit_path, 'fitResults');
    if isfield(Tfit,'fitResults'), fitResults = Tfit.fitResults; end
end

n = numel(force_dir);
Lc = nan(n,1); lk = nan(n,1); ks = nan(n,1);
Fru = nan(n,1); rload = nan(n,1);
ok  = false(n,1); min_sep = nan(n,1); max_sep = nan(n,1);
base_out = strings(n,1);
newEntries = repmat(emptyEntry(), n, 1);

% ---------------- Serial loop ----------------
for i = 1:n
    base = force_dir(i).name;

    % Skip if this curve does not have a rupture file
    if ~ismember(base, bases_with_rupture)
        continue;
    end

    % Skip if not in cache
    if ~isKey(info_map, base), continue; end
    pair = info_map(base);
    lo = pair(1); up = pair(2);
    if ~(isfinite(lo) && isfinite(up) && up > lo), continue; end

    fileName = [base, '_retract_scaled.mat_rupture.mat'];
    fullPath = fullfile(path_processed, fileName);
    if ~isfile(fullPath)
        % skip quietly; no warnings
        continue;
    end

    base_out(i) = string(fileName);

    % create offscreen axes if plotting
    if opt.SavePlots
        f = figure('Visible','off','Color','w','Position',[100 100 560 420]);
        ax = axes('Parent',f);
    else
        ax = [];
    end

    try
        % --- call your existing fitter ---
        [rload(i), Fru(i), z_model, F_model, p_fit] = fit_mFJC( ...
            path_processed, fileName, lo, up, ...
            L_cont0, T0, l_kuhn0, k_seg0, v0, ax);

        Lc(i) = p_fit(1); lk(i) = p_fit(2); ks(i) = p_fit(3);

        % --- export plot if requested ---
        if opt.SavePlots && ~isempty(ax) && isvalid(ax)
            hold(ax,'on');
            plot(ax, z_model, F_model, 'r-', 'LineWidth', 2, 'DisplayName','mFJC fit');
            title(ax, sprintf('%s | L_c=%.2g m, l_k=%.2g m, k_s=%.2g N/m', ...
                fileName, Lc(i), lk(i), ks(i)));
            xlabel(ax,'Separation (m)'); ylabel(ax,'Force (N)');
            grid(ax,'on'); legend(ax,'show');
            exportgraphics(ax, fullfile(opt.OutDir, erase(fileName,'.mat') + "_mFJC.png"), 'Resolution', 200);
            close(ancestor(ax,'figure'));
        end

        % --- store result in UI schema ---
        entry = emptyEntry();
        entry.base_str        = fileName;
        entry.initialParams   = [v0, L_cont0, T0, l_kuhn0, k_seg0];
        entry.parameters      = [Lc(i), lk(i), ks(i), T0, 1];
        entry.loading_rate    = rload(i);
        entry.unbinding_force = Fru(i);
        entry.z_model         = z_model;
        entry.F_model         = F_model;
        entry.intervals       = [lo, up];

        newEntries(i) = entry;
        ok(i) = true;
        min_sep(i) = lo; max_sep(i) = up;

    catch ME
        if ~opt.Quiet
            warning('fit_mFJC failed for %s: %s', fileName, ME.message);
        end
        if opt.SavePlots && ~isempty(ax) && isvalid(ax)
            close(ancestor(ax,'figure'));
        end
    end
end

% ---------------- Merge into fitResults ----------------
for i = 1:n
    if ~ok(i), continue; end
    e = newEntries(i);

    if ~isempty(fitResults)
        oldF = fieldnames(fitResults);
        newF = fieldnames(e);
        addToOld = setdiff(newF, oldF);
        for f = addToOld', [fitResults.(f{1})] = deal([]); end
        addToNew = setdiff(oldF, newF);
        for f = addToNew', e.(f{1}) = []; end
    end

    if opt.ReplaceExisting && ~isempty(fitResults)
        keep = ~strcmp({fitResults.base_str}, e.base_str);
        fitResults = fitResults(keep);
    end

    fitResults = [fitResults; e]; %#ok<AGROW>
end

% ---------------- Save and return ----------------
save(fit_path, 'fitResults');

results = table(base_out, ok, min_sep, max_sep, Lc, lk, ks, Fru, rload, ...
    'VariableNames', {'base_str','fit_ok','min_sep','max_sep', ...
                      'L_c','l_kuhn','k_segment','F_unbind','loading_rate'});

if ~opt.Quiet
    fprintf('Saved %d fits (%d had rupture files) to %s\n', ...
        nnz(ok), numel(bases_with_rupture), fit_path);
end
end

% -------- helper --------
function e = emptyEntry()
e = struct('base_str','', 'initialParams',[], 'parameters',[], ...
           'loading_rate',[], 'unbinding_force',[], ...
           'intervals',[], 'z_model',[], 'F_model',[]);
end
