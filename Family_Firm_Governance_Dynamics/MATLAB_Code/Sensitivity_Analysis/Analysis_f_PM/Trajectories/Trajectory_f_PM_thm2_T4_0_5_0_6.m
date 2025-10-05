%%% TRAJECTORY_F_PM_THM2_T4_0.5_0.6.M - Trajectory Analysis for Flipping Point T4 vs. F_PM (Theorem 2) %%%
clear; clc; close all;

%%% ===================================================================
%%% === 1. DEFINE SYSTEM PARAMETERS (FIXED BASELINE FOR THEOREM 2)
%%% ===================================================================
% These are the baseline parameters for the Theorem 2 scenario
phi = 0.25; m = 0.12; omega = 0.20;
alpha = 0.60; beta = 0.40;
V1 = 6; V2 = 13; V3 = 16; V4 = 10;
n = 0.05;
f1 = 3.0;
f_alpha = 1.8;

%%% ===================================================================
%%% === 2. SIMULATION CONFIGURATION
%%% ===================================================================
% The valid range for f_PM for Theorem 2 is (0, 3.09)
f_PM_values = 0.1:0.2:3.0; % Select a range of values for f_PM to test
% Define the single initial condition for this specific point
initial_condition = [0.5, 0.6];

%%% ===================================================================
%%% === 3. PLOTTING SETUP & TRAJECTORY CALCULATION
%%% ===================================================================
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 10 8]); 
hold on;
colors = parula(length(f_PM_values));
legend_handles = gobjects(length(f_PM_values), 1);

% --- Loop to plot trajectories for each 'f_PM' from the single initial point ---
for i = 1:length(f_PM_values)
    f_PM = f_PM_values(i);
    
    ode_function = @(t, y) replicator_dynamics_local(t, y, phi, m, n, omega, alpha, beta, V1, V2, V3, V4, f1, f_alpha, f_PM);
    
    [~, Y] = ode45(ode_function, [0 500], initial_condition, options);
    
    legend_label = sprintf('$f_{PM} = %.2f$', f_PM);
    
    h = plot(Y(:,1), Y(:,2), '-', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', legend_label);
    legend_handles(i) = h;
end

%%% ===================================================================
%%% === 4. PLOT ANNOTATIONS
%%% ===================================================================
% Plot the initial starting point 'T4'
plot(initial_condition(1), initial_condition(2), 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(initial_condition(1), initial_condition(2) - 0.04, '$T_4$', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

% Plot the pure-strategy equilibria with correct stability colors for Theorem 2
plot(0, 1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.1 0.4 0.8], 'MarkerEdgeColor', 'k'); % E2 is stable
text(0-0.02, 1+0.03, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.3 0], 'MarkerEdgeColor', 'k'); % E3 is stable
text(1+0.02, 0-0.03, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
color_unstable = [0.3 0.3 0.3]; 
plot(0, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', color_unstable, 'MarkerEdgeColor', 'k'); % E1 is unstable
text(0-0.02, 0-0.03, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', color_unstable, 'MarkerEdgeColor', 'k'); % E4 is unstable
text(1+0.02, 1+0.03, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

%%% ===================================================================
%%% === 5. FINAL FIGURE FORMATTING
%%% ===================================================================
axis([0 1 0 1]); axis square;
xlabel('Family Strategy, $p$ (Probability of Concentration)','Interpreter','latex', 'FontSize', 14);
ylabel('Manager Strategy, $q$ (Probability of Stewardship)','Interpreter','latex', 'FontSize', 14);
title({'Trajectory Sensitivity to $f_{PM}$ from Flipping Initial Point $T_4$ (0.5, 0.6)', '--- Theorem 2'},'Interpreter','latex', 'FontSize', 16);
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'FontName', 'Helvetica');
hold off;

% Create and format the legend
lgd = legend(legend_handles, 'Interpreter', 'latex', 'Location', 'northeastoutside');
lgd.Title.String = '\hspace{.5em}Trajectory for varying $f_{PM}$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';

% Save Figure
ax = gca;
ax.Toolbar.Visible = 'off';
outputFileName = 'Trajectory_f_PM_thm2_T4_0.5_0.6.pdf';
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);

%%% ===================================================================
%%% === LOCAL FUNCTIONS
%%% ===================================================================
function dydt = replicator_dynamics_local(~, y, phi, m, n, omega, alpha, beta, V1, V2, V3, V4, f1, f_alpha, f_PM)
    p=y(1); q=y(2); dydt=zeros(2,1);
    U_FC = q*(1-omega)*V1 + (1-q)*(1-omega-phi+m)*V2 + f1;
    U_FD = q*(alpha*(1-omega)+(1-beta)*(1-alpha)*(phi-n))*V3 + (1-q)*alpha*(1-omega-phi+m+n)*V4 + f_alpha;
    U_PS = p*omega*V1 + (1-p)*(omega*V3 + beta*(1-alpha)*(phi-n)*V3) + f_PM;
    U_PA = p*(omega+phi-m)*V2 + (1-p)*(omega+phi-m-n)*V4;
    dydt(1) = p*(1-p)*(U_FC-U_FD);
    dydt(2) = q*(1-q)*(U_PS-U_PA);
end