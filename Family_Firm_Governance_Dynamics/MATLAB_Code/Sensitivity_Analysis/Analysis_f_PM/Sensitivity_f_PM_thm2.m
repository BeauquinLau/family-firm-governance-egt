%%% SENSITIVITY_F_PM_THM2.M - SENSITIVITY ANALYSIS FOR F_PM (THEOREM 2) %%%
clear; clc; close all;

%%% ===================================================================
%%% === 1. DEFINE SYSTEM PARAMETERS (FIXED BASELINE FOR THEOREM 2)
%%% ===================================================================
% These are the baseline parameters for the Theorem 2 scenario
% NOTE: 'f_PM' is now the variable.
phi = 0.25; m = 0.12; omega = 0.20;
alpha = 0.60; beta = 0.40;
V1 = 6; V2 = 13; V3 = 16; V4 = 10;
n = 0.05;
f1 = 3.0;
f_alpha = 1.8;
f_diff = f1 - f_alpha;

%%% ===================================================================
%%% === 2. SIMULATION CONFIGURATION
%%% ===================================================================
% The valid range for f_PM was derived as (0, 3.09)
f_PM_values = 0.1:0.2:3.0; % Select a range of values for f_PM to test
p_star = zeros(size(f_PM_values));
q_star = zeros(size(f_PM_values));

%%% ===================================================================
%%% === 3. CALCULATE SADDLE POINT COORDINATES FOR EACH 'f_PM'
%%% ===================================================================
for i = 1:length(f_PM_values)
    f_PM = f_PM_values(i); % Use the loop variable for f_PM
    
    % The formulas for p* depend on f_PM
    Ap = (omega + phi - m - n) * V4 - (omega*V3 + beta*(1-alpha)*(phi-n)*V3 + f_PM);
    Bp = (omega*V1 - (omega + phi - m)*V2) - (omega*V3 + beta*(1-alpha)*(phi-n)*V3 - (omega + phi - m - n)*V4);
    p_star(i) = Ap / Bp;
    
    % The formulas for q* do not depend on f_PM, so they are constant.
    % However, to strictly follow the framework, we recalculate them in the loop.
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
colors = parula(length(f_PM_values)); 
legend_handles = gobjects(length(f_PM_values), 1);

% --- Main loop to calculate and plot manifolds for each 'f_PM' ---
for i = 1:length(f_PM_values)
    f_PM = f_PM_values(i); % Ensure f_PM is updated for the local functions
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
    
    h = plot(Y_sta_pos(:,1), Y_sta_pos(:,2), '-', 'Color', colors(i,:), 'LineWidth', 2.5, 'DisplayName', sprintf('$f_{PM}=%.2f$', f_PM));
    
    legend_handles(i) = h;
    plot(Y_sta_neg(:,1), Y_sta_neg(:,2), '-', 'Color', colors(i,:), 'LineWidth', 2.5, 'HandleVisibility','off');
    
    plot(p_val, q_val, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','k', 'MarkerSize', 8, 'HandleVisibility','off');
end

%%% ===================================================================
%%% === 5. PLOT FIXED POINTS AND FINAL FORMATTING
%%% ===================================================================
plot(0.2, 0.8, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.2, 0.8 - 0.04, '$S_7$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
plot(0.7, 0.3, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.7, 0.3 + 0.04, '$S_8$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
plot(0.5, 0.6, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.5 + 0.02, 0.6 - 0.04, '$T_4$', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex');

% Plot Fixed Points (for Theorem 2)
plot(0, 1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.1 0.4 0.8], 'MarkerEdgeColor', 'k'); % E2 is stable
text(0-0.02, 1+0.03, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.3 0], 'MarkerEdgeColor', 'k'); % E3 is stable
text(1+0.02, 0-0.03, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
color_unstable = [0.3 0.3 0.3]; 
plot(0, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', color_unstable, 'MarkerEdgeColor', 'k'); % E1 is unstable
text(0-0.02, 0-0.03, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', color_unstable, 'MarkerEdgeColor', 'k'); % E4 is unstable
text(1+0.02, 1+0.03, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% Final Figure Formatting
axis([0 1 0 1]); axis square;
xlabel('Family Strategy, $p$ (Probability of Concentration)','Interpreter','latex', 'FontSize', 14);
ylabel('Manager Strategy, $q$ (Probability of Stewardship)','Interpreter','latex', 'FontSize', 14);
title('Sensitivity Analysis for Theorem 2: Locus of the Separatrix for varying $f_{PM}$','Interpreter','latex', 'FontSize', 16);
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'FontName', 'Helvetica', 'Layer', 'top');
lgd = legend(legend_handles, 'Interpreter', 'latex', 'Location', 'northeastoutside');
lgd.Title.String = '\hspace{.5em}Separatrix for varying $f_{PM}$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';
hold off;

% Save Figure
ax = gca;
ax.Toolbar.Visible = 'off';
outputFileName = 'Sensitivity_f_PM_thm2.pdf';
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