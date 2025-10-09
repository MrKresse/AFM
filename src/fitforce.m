function [SSEc, mc, bc, SSE, m, b] = fitforce(path_incoming, base_str, varargin)
    % FITFORCE  Linear fits for constant-compliance and baseline segments.
    % Usage (auto-detect all bounds):
    %   [SSEc, mc, bc, SSE, m, b] = fitforce(path_incoming, base_str)
    %
    % Usage (override any subset; pass [] to auto-detect that bound):
    %   [SSEc, mc, bc, SSE, m, b] = fitforce(path_incoming, base_str, ...
    %       'const_compliance_min', v1, 'const_compliance_max', v2, ...
    %       'baseline_min', v3, 'baseline_max', v4, 'robust','off');
    %
    % Outputs:
    %   SSEc, mc, bc  : SSE / slope / intercept for constant-compliance fit
    %   SSE,  m,  b   : SSE / slope / intercept for baseline fit
    
    % ---------- Parse optional args ----------
    ip = inputParser;
    ip.addParameter('const_compliance_min', [], @(x) isempty(x) || isscalar(x));
    ip.addParameter('const_compliance_max', [], @(x) isempty(x) || isscalar(x));
    ip.addParameter('baseline_min',         [], @(x) isempty(x) || isscalar(x));
    ip.addParameter('baseline_max',         [], @(x) isempty(x) || isscalar(x));
    ip.addParameter('robust',            'off', @(s) ischar(s) || isstring(s)); % 'off' | 'LAR' | 'Bisquare'
    ip.parse(varargin{:});
    par = ip.Results;
    
    % ---------- Load data ----------
    approach_file = fullfile(path_incoming, base_str);
    S = load(approach_file, 'height', 'deflection');
    height     = S.height(:);
    deflection = S.deflection(:);
    
    % Clean NaNs/Infs
    good = isfinite(height) & isfinite(deflection);
    height = height(good);
    deflection = deflection(good);
    
    if numel(height) < 2
        error('fitforce:NotEnoughData', 'Need at least 2 finite data points.');
    end
    
    % ---------- Auto-detect bounds if not provided ----------
    hmin = min(height); hmax = max(height);
    
    cc_min = par.const_compliance_min; if isempty(cc_min), cc_min = hmin; end
    cc_max = par.const_compliance_max; if isempty(cc_max), cc_max = hmax; end
    bl_min = par.baseline_min;         if isempty(bl_min), bl_min = hmin; end
    bl_max = par.baseline_max;         if isempty(bl_max), bl_max = hmax; end
    
    % Ensure ordering (swap if user inverted)
    if cc_min > cc_max, [cc_min, cc_max] = deal(cc_max, cc_min); end
    if bl_min > bl_max, [bl_min, bl_max] = deal(bl_max, bl_min); end
    
    % ---------- Build indices ----------
    idx_compliance = (height >= cc_min) & (height <= cc_max);
    idx_baseline   = (height >= bl_min) & (height <= bl_max);
    
    if nnz(idx_compliance) < 2
        error('fitforce:ComplianceWindow', 'Constant-compliance window has <2 points.');
    end
    if nnz(idx_baseline) < 2
        error('fitforce:BaselineWindow', 'Baseline window has <2 points.');
    end
    
    xC = height(idx_compliance);
    yC = deflection(idx_compliance);
    xB = height(idx_baseline);
    yB = deflection(idx_baseline);
    
    % ---------- Fit (poly1) ----------
    % Robust option uses Curve Fitting Toolbox robust fitting (optional)
    fitArgs = {};
    if ~strcmpi(par.robust, 'off')
        fitArgs = {'Robust', par.robust}; % 'LAR' or 'Bisquare'
    end
    
    [fitC, gofC] = fit(xC, yC, 'poly1', fitArgs{:});
    SSEc = gofC.sse;
    mc   = fitC.p1;
    bc   = fitC.p2;
    
    [fitB, gofB] = fit(xB, yB, 'poly1', fitArgs{:});
    SSE = gofB.sse;
    m   = fitB.p1;
    b   = fitB.p2;
end
