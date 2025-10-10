function [loading_rate, unbinding_force, z_model, F_model, p_fit] = fit_mFJC(path_processed, base_str, min_sep, max_sep, l_cont, T, l_kuhn, k_segment, v, ax)
    % Compute thermal energy kBT in J
    kBT = 1.380649e-23 * T;

    % Load the data
    data = load(fullfile(path_processed, base_str));
    F = data.scaled_deflection(:); % Force data (N)
    z = data.shifted_height(:);    % Extension data (m)

    % Select the fitting interval
    inds  = (z <= max_sep) & (z >= min_sep) & isfinite(z) & isfinite(F);
    F_fit = F(inds);
    z_fit = z(inds);

    % Define the modified FJC model:
    % p(1): L_c, p(2): l_k, p(3): k_segment
    mFJC = @(p, F_val) p(1) .* (coth(F_val .* p(2) ./ kBT) - kBT./(F_val .* p(2))) ...
                       .* (1 + F_val./(p(3).*p(2)));

    % Initial guess and bounds
    p0 = [l_cont, l_kuhn, k_segment];
    lb = [0, 0, 0];
    ub = [Inf, Inf, Inf];

    % Nonlinear least-squares fit (fit in positive force magnitude space)
    options = optimoptions('lsqcurvefit', 'Display','iter', ...
        'FunctionTolerance',1e-30, 'OptimalityTolerance',1e-25, 'StepTolerance',1e-25);
    [p_fit, ~, ~, ~, ~] = lsqcurvefit(mFJC, p0, -F_fit, z_fit, lb, ub, options);

    % Build model curve (use positive forces for computation)
    F_max_mag = max(1e-12, -min(F));  % ensure > 0, avoids empty/zero ranges
    F_range   = linspace(0, F_max_mag, 100);
    z_model   = mFJC(p_fit, F_range);
    F_model   = -F_range;             % revert to original sign convention

    % Optional plotting (only if a valid axes handle is provided)
    doPlot = (nargin >= 10) && ~isempty(ax) && isgraphics(ax,'axes');
    if doPlot
        cla(ax);
        plot(ax, z, F, 'b', 'DisplayName', 'Retract'); hold(ax, 'on');
        plot(ax, z_model, F_model, 'r-', 'LineWidth', 2, 'DisplayName', 'mFJC model');
        xlabel(ax, 'Separation (m)'); ylabel(ax, 'Force (N)');
        legend(ax, 'show'); grid(ax, 'on'); hold(ax, 'off');
    end

    % Unbinding force from experimental fit interval
    [unbinding_force, ~] = min(F_fit);

    % Effective stiffness at the rupture force (dF/dz from model)
    dF_dz_model = gradient(F_model, z_model);
    [~, i_model] = min(abs(F_model - unbinding_force));
    k_eff = dF_dz_model(i_model);

    % Loading rate
    loading_rate = v * k_eff;  % N/s
end
