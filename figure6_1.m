clear; clc; close all

%% ===== PARAMETERS =====
A = 0.92;      % retention
b = 0.02;      % error sensitivity
e = 45;        % error clamp (deg)
n_cycles = 40; % learning cycles

%% ===== ENGAGEMENT PROFILES =====
i = 0:n_cycles-1;

% ST
m_ST = ones(1, n_cycles);

% DT Option A: Engagement boost that decays to baseline
delta = 0.35;
rho = 0.9;
m_DT_decay = 1 + delta * rho.^i;

% DT Option B: Partial boost in engagement
m_DT_partial = ones(1, n_cycles) * (1 + 0.10);

% DTF: Sustained boost in engagement
m_DTF = ones(1, n_cycles) * (1 + delta);

%% ===== RUN SIMULATIONS =====
x_ST         = sim_adapt(A, b, e, m_ST);
x_DT_decay   = sim_adapt(A, b, e, m_DT_decay);
x_DT_partial = sim_adapt(A, b, e, m_DT_partial);
x_DTF        = sim_adapt(A, b, e, m_DTF);

%% ===== FIGURE 6-1 =====
t = 1:n_cycles;
group_colors = struct('ST',hex2rgb("#3182bd"), 'DT', hex2rgb("#c51b8a"), 'DTF', hex2rgb("#ff8c12"));

f = figure("Color", "w", "Units", "inches", "OuterPosition", [3 2 6.5 3.25]);
theme(f, "light");

% Panel a: Hand angle
subplot(1,2,1); hold on
plot(t+10, x_ST,         'Color', group_colors.ST,  'LineWidth', 2);
plot(t+10, x_DT_decay,   'Color', group_colors.DT,  'LineWidth', 2);
plot(t+10, x_DT_partial, ':',     'Color', group_colors.DT,  'LineWidth', 1.5);
plot(t+10, x_DTF,        'Color', group_colors.DTF, 'LineWidth', 2);
xlim([10, 51]);
ylim([0, 15]);
xlabel('cycle');
ylabel('hand angle (°)')
legend({'ST','DT (decay)','DT (partial)','DT_{F} (sustained)'}, 'Location', 'southeast', "Box", "off");
annotation("textbox", [.03, .97, .0, 0], "String", "a", "FontName", "Arial", "FontSize", 14, "FontWeight", "bold", "LineStyle", "none");

% Panel b: Change in hand angle from ST
subplot(1,2,2); hold on;
plot(t+10, x_DT_decay   - x_ST, '-',  'Color', group_colors.DT,  'LineWidth', 1.5)
plot(t+10, x_DT_partial - x_ST, ':',  'Color', group_colors.DT,  'LineWidth', 1.5)
plot(t+10, x_DTF        - x_ST, '-',  'Color', group_colors.DTF, 'LineWidth', 1.5)
yline(0, 'k:', 'LineWidth', 0.5);
xlim([10, 51]);
ylim([0, 6]);
xlabel('cycle');
ylabel('\Delta hand angle (°)')
legend({'DT (decay)','DT (partial)','DT_{F} (sustained)'}, 'Location', 'northeast', "Box", "off");
annotation("textbox", [.48, .97, .0, 0], "String", "b", "FontName", "Arial", "FontSize", 14, "FontWeight", "bold", "LineStyle", "none");

% Save figure
exportgraphics(f, fullfile("figures", "figure6_1.tif"), 'Resolution', 600);

%% Supporting functions
% Simulate single-state model: x_{n+1} = A*x_n + b*m(n)*e
function x = sim_adapt(A, b, e, m)
    n = numel(m);
    x = zeros(1, n);
    for k = 1:n-1
        x(k+1) = A*x(k) + b*m(k)*e;
    end
end
