function [SSEc, mc, bc, SSE, m, b]=fitforce(path_incoming,base_str, const_compliance_max,baseline_min)
       
    % Define the filenames based on the base_str
    approach_file = fullfile(path_incoming, base_str);
     % Load the approach and deflection data from the saved .mat files
    load(approach_file, 'height', 'deflection');

    const_compliance_min=min(height);
    baseline_max=max(height);
    % Fit a linear function to the constant compliance interval
    idx_compliance = (height >= const_compliance_min) & (height <= const_compliance_max);
    x_compliance = height(idx_compliance);
    y_compliance = deflection(idx_compliance);
    [fit_compliance, gof_compliance]= fit(x_compliance, y_compliance, 'poly1');  % Fit a linear function
    SSEc = gof_compliance.sse;
    mc = fit_compliance.p1;
    bc = fit_compliance.p2;

    % Fit a linear function to the baseline interval
    idx_baseline = (height >= baseline_min) & (height <= baseline_max);
    x_baseline = height(idx_baseline);
    y_baseline = deflection(idx_baseline);
    [fit_baseline, gof_baseline] = fit(x_baseline, y_baseline, 'poly1');  % Fit a linear function
    SSE = gof_baseline.sse;
    m = fit_baseline.p1;
    b = fit_baseline.p2;
end