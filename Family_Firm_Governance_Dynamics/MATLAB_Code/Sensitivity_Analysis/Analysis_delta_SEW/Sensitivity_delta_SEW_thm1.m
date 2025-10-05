%%% SENSITIVITY_DELTA_SEW_THM1.M - SENSITIVITY ANALYSIS FOR PARAMETER DELTA_SEW %%%
clear; clc; close all;

%%% ===================================================================
%%% === 1. DEFINE SYSTEM PARAMETERS (FIXED BASELINE)
%%% ===================================================================
% These are the baseline parameters for the Theorem 1 scenario
phi = 0.3; m = 0.04; omega = 0.2; alpha = 0.5; beta = 0.5;
V1 = 10; V2 = 5; V3 = 8; V4 = 15;
f_PM = 1.2;
n = 0.16;
f_alpha = 1.0; % We fix f_alpha and vary f1 to change the SEW gap

%%% ===================================================================
%%% === 2. SIMULATION CONFIGURATION
%%% ===================================================================
% The valid range for delta_SEW was derived as [0, 2.55)
delta_SEW_values = 0.1:0.2:2.5; % Select a range of values for delta_SEW to test
p_star = zeros(size(delta_SEW_values));
q_star = zeros(size(delta_SEW_values));

%%% ===================================================================
%%% === 3. CALCULATE SADDLE POINT COORDINATES FOR EACH 'delta_SEW'
%%% ===================================================================
for i = 1:length(delta_SEW_values)
    delta_SEW = delta_SEW_values(i);
    f1 = f_alpha + delta_SEW; % Calculate f1 based on the current delta_SEW
    f_diff = delta_SEW; % Use the loop variable directly
    
    % The formulas for p* do not depend on f(1) or f(alpha), so Ap and Bp are constant.
    % However, to strictly follow the framework, we recalculate them in the loop.
    Ap = (omega + phi - m - n) * V4 - (omega*V3 + beta*(1-alpha)*(phi-n)*V3 + f_PM);
    Bp = (omega*V1 - (omega + phi - m)*V2) - (omega*V3 + beta*(1-alpha)*(phi-n)*V3 - (omega + phi - m - n)*V4);
    p_star(i) = Ap / Bp;
    
    % The formulas for q* depend on f_diff (our delta_SEW)
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
colors = parula(length(delta_SEW_values)); 
legend_handles = gobjects(length(delta_SEW_values), 1);

% --- Main loop to calculate and plot manifolds for each 'delta_SEW' ---
for i = 1:length(delta_SEW_values)
    delta_SEW = delta_SEW_values(i);
    f1 = f_alpha + delta_SEW;
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
    
    % Updated DisplayName for the legend
    h = plot(Y_sta_pos(:,1), Y_sta_pos(:,2), '-', 'Color', colors(i,:), 'LineWidth', 2.5, 'DisplayName', sprintf('$\\Delta_{SEW}=%.2f$', delta_SEW));
    
    legend_handles(i) = h;
    plot(Y_sta_neg(:,1), Y_sta_neg(:,2), '-', 'Color', colors(i,:), 'LineWidth', 2.5, 'HandleVisibility','off');
    
    plot(p_val, q_val, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','k', 'MarkerSize', 8, 'HandleVisibility','off');
end

%%% ===================================================================
%%% === 5. PLOT FIXED POINTS AND FINAL FORMATTING
%%% ===================================================================
% Plotting the five initial points
plot(0.1, 0.1, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.1, 0.1 - 0.04, '$S_1$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
plot(0.6, 0.5, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.6, 0.5 + 0.04, '$S_2$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
plot(0.2, 0.4, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(0.2 + 0.02, 0.4 - 0.04, '$T_1$', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex');

% Plot Fixed Points
plot(0, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.3 0], 'MarkerEdgeColor', 'k');
text(0-0.02, 0-0.03, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.1 0.4 0.8], 'MarkerEdgeColor', 'k');
text(1+0.02, 1+0.03, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
color_unstable = [0.3 0.3 0.3]; 
plot(0, 1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', color_unstable, 'MarkerEdgeColor', 'k');
text(0-0.02, 1+0.03, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', color_unstable, 'MarkerEdgeColor', 'k');
text(1+0.02, 0-0.03, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% Final Figure Formatting
axis([0 1 0 1]); axis square;
xlabel('Family Strategy, $p$ (Probability of Concentration)','Interpreter','latex', 'FontSize', 14);
ylabel('Manager Strategy, $q$ (Probability of Stewardship)','Interpreter','latex', 'FontSize', 14);
title('Sensitivity Analysis for Theorem 1: Locus of the Separatrix for varying $\Delta_{SEW}$','Interpreter','latex', 'FontSize', 16);
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'FontName', 'Helvetica', 'Layer', 'top');
lgd = legend(legend_handles, 'Interpreter', 'latex', 'Location', 'northeastoutside');
lgd.Title.String = '\hspace{.5em}Separatrix for varying $\Delta_{SEW}$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';
hold off;

% Save Figure
ax = gca;
ax.Toolbar.Visible = 'off';
outputFileName = 'Sensitivity_delta_SEW_thm1.pdf';
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