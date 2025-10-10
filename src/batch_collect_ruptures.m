function out = batch_collect_ruptures(force_dir, path_processed, lower, upper, varargin)
% BATCH_COLLECT_RUPTURES  Identify and rank rupture events across all curves (no plotting).
%
% Required inputs:
%   force_dir      : struct array (from dir(...)) with .name for each curve
%   path_processed : folder containing '<base>_retract_scaled.mat' or
%                    '<base>_retract_scaled.mat_rupture.mat' if 'use_rupture_files' = true
%   lower, upper   : preferred separation window [m]
%
% Optional name-value arguments (defaults in []):
%   'dz_bin'          [1e-9]   bin width in separation (m)
%   'min_pts_per_bin' [1]     minimum raw points per bin
%   'use_movmedian'   [5]     odd window for movmedian (0 disables)
%   'strict_window'   [true]  require z_pre>=lower & z_post<=upper
%   'mask_zmax'       [1e-7]  base mask: only use z <= mask_zmax
%   'use_rupture_files' [false] load '<base>_retract_scaled.mat_rupture.mat' instead of default
%   'write_cache'       [true]  save rupture_info.mat cache next to processed data
%   'cache_file'        [fullfile(path_processed,'rupture_info.mat')]
%
% Output struct 'out' (all SI units: N, m):
%   .results(:) : per-curve struct with fields:
%       curve_id, base_str, found, z_pre, z_post, F_pre, F_post, dF, n_points, n_bins
%   .dF_keep, .z_keep, .idx_curve
%   .sorted_idx : curve indices sorted by descending ΔF
%   .sorted_dF  : sorted ΔF values (N)
%   .table      : MATLAB table view of results
%
% Also (if write_cache): writes 'rupture_info.mat' with per-curve metadata:
%   rupture_info(i) = struct('base_str','...','found',true,'z_pre',...,'z0_left',...,
%                            'z_post',...,'F_pre',...,'F_post',...,'dF',...,
%                            'n_points',...,'n_bins',...,'params',p,'timestamp',datetime)

% ---------- Parse inputs ----------
ip = inputParser;
ip.addParameter('dz_bin',          1e-9);
ip.addParameter('min_pts_per_bin', 1);
ip.addParameter('use_movmedian',   5);
ip.addParameter('strict_window',   true);
ip.addParameter('mask_zmax',       1e-7);
ip.addParameter('use_rupture_files', false);
ip.addParameter('write_cache',     true);
ip.addParameter('cache_file',      fullfile(path_processed,'rupture_info.mat'));
ip.parse(varargin{:});
p = ip.Results;

num_curves = numel(force_dir);

% ---------- Preallocate ----------
results(num_curves,1) = struct('curve_id',NaN,'base_str','','found',false,...
    'z_pre',NaN,'z_post',NaN,'F_pre',NaN,'F_post',NaN,'dF',NaN,...
    'n_points',0,'n_bins',0);

dF_keep   = nan(num_curves,1);
z_keep    = nan(num_curves,1);
idx_curve = false(num_curves,1);

rupture_info = struct('base_str',{},'found',{},'z_pre',{},'z0_left',{},'z_post',{},...
                      'F_pre',{},'F_post',{},'dF',{},'n_points',{},'n_bins',{},...
                      'params',{},'timestamp',{});

% ---------- Batch process each curve ----------
for cid = 1:num_curves
    base_str = force_dir(cid).name;

    % Choose file variant
    if p.use_rupture_files
        file = fullfile(path_processed, [base_str, '_retract_scaled.mat_rupture.mat']);
    else
        file = fullfile(path_processed, [base_str, '_retract_scaled.mat']);
    end
    if ~isfile(file), continue; end

    S = load(file);
    if ~isfield(S,'shifted_height') || ~isfield(S,'scaled_deflection'), continue; end

    z = S.shifted_height(:);
    F = S.scaled_deflection(:);
    mask = (z <= p.mask_zmax) & isfinite(z) & isfinite(F);
    if ~any(mask), continue; end

    z_use = z(mask);  F_use = F(mask);
    [z_sorted, ord] = sort(z_use);
    F_sorted = F_use(ord);
    if p.use_movmedian > 0
        F_sorted = movmedian(F_sorted, p.use_movmedian);
    end

    % ---- Bin by separation (medians) ----
    zmin = max(min(z_sorted), lower);
    zmax = min(max(z_sorted), upper);
    if ~(zmax > zmin), continue; end

    edges = zmin:p.dz_bin:zmax;
    if edges(end) < zmax, edges(end+1) = zmax; end
    b = discretize(z_sorted, edges);
    if all(isnan(b)), continue; end
    nb = max(b(~isnan(b)));

    F_bin = nan(nb,1); z_bin = nan(nb,1); nbin = zeros(nb,1);
    for k = 1:nb
        m = (b == k);
        if any(m)
            F_bin(k) = median(F_sorted(m));
            z_bin(k) = median(z_sorted(m));
            nbin(k)  = sum(m);
        end
    end
    valid = isfinite(F_bin) & isfinite(z_bin) & (nbin >= p.min_pts_per_bin);
    F_bin = F_bin(valid); z_bin = z_bin(valid);
    if numel(F_bin) < 2, continue; end

    % ---- Compute ΔF between bins ----
    dF = diff(F_bin);
    zL = z_bin(1:end-1); zR = z_bin(2:end);
    FL = F_bin(1:end-1); FR = F_bin(2:end);

    pos_jump = (dF > 0);
    neg_pre  = (FL < 0);

    if p.strict_window
        in_win = (zL >= lower) & (zR <= upper);
    else
        in_win = true(size(zL));
    end

    cand = find(pos_jump & neg_pre & in_win);
    if isempty(cand)
        pos_all = find(dF > 0);
        if isempty(pos_all), continue; end
        [~, imax] = max(dF(pos_all));
        k = pos_all(imax);
    else
        [~, imax] = max(dF(cand));
        k = cand(imax);
    end

    % ---- Record result ----
    results(cid).curve_id = cid;
    results(cid).base_str = base_str;
    results(cid).found    = true;
    results(cid).z_pre    = zL(k);
    results(cid).z_post   = zR(k);
    results(cid).F_pre    = FL(k);
    results(cid).F_post   = FR(k);
    results(cid).dF       = dF(k);
    results(cid).n_points = numel(z_use);
    results(cid).n_bins   = numel(F_bin);

    dF_keep(cid) = dF(k);
    z_keep(cid)  = zL(k);
    idx_curve(cid) = true;

    % ---- Compute first zero-crossing LEFT of z_pre ----
    z0_left = local_first_zero_left(z_sorted, F_sorted, zL(k), p.dz_bin);

    % ---- Append to cache ----
    rupture_info(end+1,1) = struct( ...
        'base_str', base_str, ...
        'found', true, ...
        'z_pre', zL(k), ...
        'z0_left', z0_left, ...
        'z_post', zR(k), ...
        'F_pre', FL(k), ...
        'F_post', FR(k), ...
        'dF', dF(k), ...
        'n_points', numel(z_use), ...
        'n_bins', numel(F_bin), ...
        'params', p, ...
        'timestamp', datetime('now'));
end

% ---------- Sort all by jump magnitude ----------
valid_idx = find(idx_curve);
if ~isempty(valid_idx)
    [sorted_dF, order] = sort(dF_keep(valid_idx), 'descend');
    sorted_idx = valid_idx(order);
else
    sorted_dF = [];
    sorted_idx = [];
end

% ---------- Build table & output ----------
T = struct2table(results);

out = struct( ...
    'results',     results, ...
    'idx_curve',   idx_curve, ...
    'dF_keep',     dF_keep, ...
    'z_keep',      z_keep, ...
    'sorted_idx',  sorted_idx, ...
    'sorted_dF',   sorted_dF, ...
    'table',       T );

%fprintf('Largest ΔF = %.2f pN; Median = %.2f pN.\n', ...
%        (isempty(sorted_dF) * NaN) + (~isempty(sorted_dF) * max(sorted_dF))*1e12, ...
%        (isempty(sorted_dF) * NaN) + (~isempty(sorted_dF) * median(sorted_dF))*1e12 );

% ---------- Write cache ----------
if p.write_cache
    try
        save(p.cache_file, 'rupture_info');
    catch ME
        warning('Could not save rupture_info cache: %s', ME.message);
    end
end

end % main function

% ===== helper =====
function z0 = local_first_zero_left(z_sorted, F_sorted, z_ref, dz_hint)
% Find first zero-crossing to the LEFT of z_ref using sorted data
z0 = NaN;
left = z_sorted < z_ref;
if nnz(left) < 2
    z0 = min(z_sorted);
    return;
end
zl = z_sorted(left); Fl = F_sorted(left);
sgn = sign(Fl);
idx = find(diff(sgn) ~= 0);
if ~isempty(idx)
    c = idx(end);
    z1 = zl(c); z2 = zl(c+1); f1 = Fl(c); f2 = Fl(c+1);
    if f2 ~= f1
        z0 = z1 + (0 - f1) * (z2 - z1) / (f2 - f1);
    end
end
if ~isfinite(z0), z0 = max(min(z_sorted), z_ref - 5*dz_hint); end
if z0 >= z_ref,   z0 = z_ref - max(5*dz_hint, eps(z_ref));   end
end
