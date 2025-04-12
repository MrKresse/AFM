% mfjc_parameter_study.m
% This script investigates the behavior of the mFJC model by varying individual parameters.
% Define base parameters
L_cont_base = 1e-8;      % Base contour length (m)
l_kuhn_base = 5e-10;     % Base Kuhn length (m)
k_seg_base = 2;       % Base segment elasticity (N/m)
T = 300;                 % Temperature (K)
kBT = 1.380649e-23 * T;  % Thermal energy (J)

% Define force range for the model (use positive force values for computation)
F_range = linspace(1e-12, 1e-9, 200); % Avoid F=0 to prevent division by zero

% Define the mFJC model as an anonymous function
% p(1): contour length L_c, p(2): Kuhn length l_k, p(3): segment elasticity k_seg.
mFJC = @(p, F_val) p(1) .* (coth(F_val .* p(2) ./ kBT) - kBT./(F_val .* p(2))) .* (1 + F_val./(p(3).*p(2)));

%% Plot baseline model output
p_base = [L_cont_base, l_kuhn_base, k_seg_base];
z_model_base = mFJC(p_base, F_range);
figure;
plot(z_model_base, -F_range, 'k-', 'LineWidth', 2);
hold on;
title('mFJC Model Behavior');
xlabel('Extension (m)');
ylabel('Force (N)');
legendEntries = {'Baseline'};

%% Vary Contour Length (L_cont) shifts the turning point to the right
L_cont_vals = [0.5, 1, 2] * L_cont_base;
colors = {'r', 'g', 'b'};
for i = 1:length(L_cont_vals)
    p = [L_cont_vals(i), l_kuhn_base, k_seg_base];
    z_model = mFJC(p, F_range);
plot(z_model, -F_range, [colors{i} '--'], 'LineWidth', 1.5);

    legendEntries{end+1} = sprintf('L_{cont} = %.2e m', L_cont_vals(i));
end

%% Vary Kuhn Length (l_kuhn) decreases the skewness
l_kuhn_vals = [0.5, 1, 2] * l_kuhn_base;
for i = 1:length(l_kuhn_vals)
    p = [L_cont_base, l_kuhn_vals(i), k_seg_base];
    z_model = mFJC(p, F_range);
 plot(z_model, -F_range, [colors{i} '--'], 'LineWidth', 1.5);

    legendEntries{end+1} = sprintf('l_{kuhn} = %.2e m', l_kuhn_vals(i));
end

%% Vary Segment Elasticity (k_seg) increases slope after turning point
k_seg_vals = [0.5, 1, 2] * k_seg_base;
for i = 1:length(k_seg_vals)
    p = [L_cont_base, l_kuhn_base, k_seg_vals(i)];
    z_model = mFJC(p, F_range);
plot(z_model, -F_range, [colors{i} '--'], 'LineWidth', 1.5);

    legendEntries{end+1} = sprintf('k_{seg} = %.2e N/m', k_seg_vals(i));
end

legend(legendEntries, 'Location', 'best');
grid on;
hold off;
