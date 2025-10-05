%%% SENSITIVITY_N_THM2.M - SENSITIVITY ANALYSIS FOR PARAMETER N (THEOREM 2) %%%
clear; clc; close all;

%%% ===================================================================
%%% === 1. DEFINE SYSTEM PARAMETERS (FIXED BASELINE FOR THEOREM 2)
%%% ===================================================================
% These are the baseline parameters for the Theorem 2 scenario
phi = 0.25; m = 0.12; omega = 0.20;
alpha = 0.60; beta = 0.40;
V1 = 6; V2 = 13; V3 = 16; V4 = 10;
f1 = 3.0; f_alpha = 1.8; f_PM = 1.8;
f_diff = f1 - f_alpha; % Pre-calculate for efficiency

%%% ===================================================================
%%% === 2. SIMULATION CONFIGURATION
%%% ===================================================================
% The valid range for n was derived as (0, 0.13)
n_values = 0.01:0.01:0.12; % Select a range of values for n to test
p_star = zeros(size(n_values));
q_star = zeros(size(n_values));

%%% ===================================================================
%%% === 3. CALCULATE SADDLE POINT COORDINATES FOR EACH 'n'
%%% ===================================================================
for i = 1:length(n_values)
    n = n_values(i);
    
    Ap = (omega + phi - m - n) * V4 - (omega*V3 + beta*(1-alpha)*(phi-n)*V3 + f_PM);
    Bp = (omega*V1 - (omega + phi - m)*V2) - (omega*V3 + beta*(1-alpha)*(phi-n)*V3 - (omega + phi - m - n)*V4);
    p_star(i) = Ap / Bp;
    Aq = alpha * (1 - omega - phi + m + n) * V4 - (1 - omega - phi + m) * V2 - f_diff;
    Bq = (1 - omega) * V1 - (1 - omega - phi + m) * V2 - (alpha*(1-omega) + (1-beta)*(1-alpha)*(phi-n))*V3 + alpha*(1-omega-phi+m+n)*V4;
    q_star(i) = Aq / Bq;
end

%%% ===================================================================
%%% === 4. PLOTTING SETUP & MANIFOLD CALCULATION
%%% ===================================================================
eps = 1e-5;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 10 8]); 
hold on;
colors = parula(length(n_values)); 
legend_handles = gobjects(length(n_values), 1);
% --- Main loop to calculate and plot manifolds for each 'n' ---
for i = 1:length(n_values)
    n = n_values(i);
    p_val = p_star(i);
    q_val = q_star(i);
    
    ode_func_handle = @(t, y) replicator_dynamics_local(t, y, phi, m, n, omega, alpha, beta, V1, V2, V3, V4, f1, f_alpha, f_PM);
    
    J_E5 = calculate_jacobian_local(p_val, q_val, phi, m, n, omega, alpha, beta, V1, V2, V3, V4, f1, f_alpha, f_PM);
    [V_eig, D_eig] = eig(J_E5);
    
    [~, idx_neg] = min(real(diag(D_eig)));
    v_stable = real(V_eig(:, idx_neg));
    init_stable_pos = [p_val; q_val] + eps * v_stable;
    init_stable_neg = [p_val; q_val] - eps * v_stable;
    [~, Y_sta_pos] = ode45(ode_func_handle, [0 -500], init_stable_pos, options);
    [~, Y_sta_neg] = ode45(ode_func_handle, [0 -500], init_stable_neg, options);
    
    h = plot(Y_sta_pos(:,1), Y_sta_pos(:,2), '-', 'Color', colors(i,:), 'LineWidth', 2.5, 'DisplayName', sprintf('$n=%.2f$', n));
    
    legend_handles(i) = h;
    plot(Y_sta_neg(:,1), Y_sta_neg(:,2), '-', 'Color', colors(i,:), 'LineWidth', 2.5, 'HandleVisibility','off');
    
    plot(p_val, q_val, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','k', 'MarkerSize', 8, 'HandleVisibility','off');
end

%%% ===================================================================
%%% === 5. PLOT FIXED POINTS AND FINAL FORMATTING
%%% ===================================================================

% Plot the five specified initial points with a red square marker
plot(0.9237, 0.8840, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.9237 - 0.04, 0.8840, '$R_2$', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex');

plot(0.6, 0.4, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.6, 0.4 - 0.04, '$S_4$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

plot(0.4, 0.7, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.4 + 0.02, 0.7 - 0.04, '$S_3$', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex');

plot(0.3, 0.4, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.3 + 0.02, 0.4 - 0.04, '$T_3$', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex');

plot(0.975740, 0.942642, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.975740 - 0.04, 0.942642, '$T_4$', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex');

% Plot Fixed Points
% E2 and E3 are now the stable points (ESS)
plot(0, 1, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.1 0.4 0.8], 'MarkerEdgeColor', 'k'); % Blue for E2
text(0-0.02, 1+0.03, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.8 0.3 0], 'MarkerEdgeColor', 'k'); % Orange for E3
text(1+0.02, 0-0.03, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
% E1 and E4 are now the unstable points
color_unstable = [0.3 0.3 0.3]; 
plot(0, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', color_unstable, 'MarkerEdgeColor', 'k');
text(0-0.02, 0-0.03, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', color_unstable, 'MarkerEdgeColor', 'k');
text(1+0.02, 1+0.03, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% Final Figure Formatting
axis([0 1 0 1]); axis square;
xlabel('Family Strategy, $p$ (Probability of Concentration)','Interpreter','latex', 'FontSize', 14);
ylabel('Manager Strategy, $q$ (Probability of Stewardship)','Interpreter','latex', 'FontSize', 14);
title('Sensitivity Analysis for Theorem 2: Locus of the Separatrix for varying $n$','Interpreter','latex', 'FontSize', 16);
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'FontName', 'Helvetica', 'Layer', 'top');
lgd = legend(legend_handles, 'Interpreter', 'latex', 'Location', 'northeastoutside');
lgd.Title.String = '\hspace{.5em}Separatrix for varying $n$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';
hold off;

% Save Figure
ax = gca;
ax.Toolbar.Visible = 'off';
outputFileName = 'Sensitivity_n_thm2.pdf';
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

function J = calculate_jacobian_local(p, q, phi, m, n, omega, alpha, beta, V1, V2, V3, V4, f1, f_alpha, f_PM)
    U_FC = q*(1-omega)*V1 + (1-q)*(1-omega-phi+m)*V2 + f1;
    U_FD = q*(alpha*(1-omega)+(1-beta)*(1-alpha)*(phi-n))*V3 + (1-q)*alpha*(1-omega-phi+m+n)*V4 + f_alpha;
    U_PS = p*omega*V1 + (1-p)*(omega*V3 + beta*(1-alpha)*(phi-n)*V3) + f_PM;
    U_PA = p*(omega+phi-m)*V2 + (1-p)*(omega+phi-m-n)*V4;
    U_F_diff = U_FC - U_FD;
    U_P_diff = U_PS - U_PA;
    dUF_dq = ((1-omega)*V1-(1-omega-phi+m)*V2) - ((alpha*(1-omega)+(1-beta)*(1-alpha)*(phi-n))*V3 - alpha*(1-omega-phi+m+n)*V4);
    dUP_dp = (omega*V1 - (omega*V3 + beta*(1-alpha)*(phi-n)*V3)) - ((omega+phi-m)*V2 - (omega+phi-m-n)*V4);
    J11 = (1-2*p)*U_F_diff;
    J12 = p*(1-p)*dUF_dq;
    J21 = q*(1-q)*dUP_dp;
    J22 = (1-2*q)*U_P_diff;
    J = [J11, J12; J21, J22];
end