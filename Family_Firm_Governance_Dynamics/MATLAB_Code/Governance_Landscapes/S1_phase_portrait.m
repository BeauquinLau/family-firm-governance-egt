% S1_phase_portrait.m - Publication-Quality Numerical Simulation for Theorem 1
clear; clc; close all;

%%% ===================================================================
%%% === 1. MODEL PARAMETERS
%%% ===================================================================
phi = 0.3; m = 0.04; n = 0.16; omega = 0.2; alpha = 0.5; beta = 0.5;
V1 = 10; V2 = 5; V3 = 8; V4 = 15; f1 = 2.5; f_alpha = 1.0; f_PM = 1.2;

%%% ===================================================================
%%% === 2. SADDLE POINT (E5) CALCULATION
%%% ===================================================================
Ap = (omega + phi - m - n) * V4 - (omega*V3 + beta*(1-alpha)*(phi-n)*V3 + f_PM);
Bp = (omega*V1 - (omega + phi - m)*V2) - (omega*V3 + beta*(1-alpha)*(phi-n)*V3 - (omega + phi - m - n)*V4);
p_star = Ap / Bp;
f_diff = f1 - f_alpha;
Aq = alpha * (1 - omega - phi + m + n) * V4 - (1 - omega - phi + m) * V2 - f_diff;
Bq = (1 - omega) * V1 - (1 - omega - phi + m) * V2 - (alpha*(1-omega) + (1-beta)*(1-alpha)*(phi-n))*V3 + alpha*(1-omega-phi+m+n)*V4;
q_star = Aq / Bq;
fprintf('Saddle point E5 is at (p*, q*) = (%.4f, %.4f)\n', p_star, q_star);

%%% ===================================================================
%%% === 3. SIMULATION SETUP
%%% ===================================================================
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 8 8]); 
hold on;
ode_function = @(t, y) replicator_dynamics(t, y, phi, m, n, omega, alpha, beta, V1, V2, V3, V4, f1, f_alpha, f_PM);

%%% ===================================================================
%%% === 4. SEPARATRIX CALCULATION & DATA PURIFICATION
%%% ===================================================================
eps = 1e-5; 
J_E5 = calculate_jacobian(p_star, q_star, phi, m, n, omega, alpha, beta, V1, V2, V3, V4, f1, f_alpha, f_PM);
[V_eig, D_eig] = eig(J_E5);
lambda = diag(D_eig);
[~, idx_neg] = min(real(lambda));
v_stable = real(V_eig(:, idx_neg));
init_stable_pos = [p_star; q_star] + eps * v_stable;
init_stable_neg = [p_star; q_star] - eps * v_stable;
[~, Y_sep_pos] = ode45(ode_function, [0 -1000], init_stable_pos, options);
[~, Y_sep_neg] = ode45(ode_function, [0 -1000], init_stable_neg, options);
separatrix_pts = [Y_sep_pos(end:-1:1, :); [p_star, q_star]; Y_sep_neg];
% --- Data Purification Pipeline ---
valid_indices = all(isfinite(separatrix_pts), 2);
separatrix_pts = separatrix_pts(valid_indices, :);
[~, unique_indices] = uniquetol(separatrix_pts, 1e-6, 'ByRows', true);
separatrix_pts = separatrix_pts(sort(unique_indices), :);
valid_indices = all(separatrix_pts >= 0, 2) & all(separatrix_pts <= 1, 2);
separatrix_pts = separatrix_pts(valid_indices, :);

%%% ===================================================================
%%% === 5. PLOT SHADED BASINS OF ATTRACTION
%%% ===================================================================
color_E4_basin = [0.85, 0.9, 0.98]; 
color_E1_basin = [0.98, 0.92, 0.82];
if ~isempty(separatrix_pts)
    vertx_E4 = [separatrix_pts(:,1); 1; 0; 0];
    verty_E4 = [separatrix_pts(:,2); 1; 1; 1];
    patch(vertx_E4, verty_E4, color_E4_basin, 'EdgeColor', 'none');
    vertx_E1 = [separatrix_pts(:,1); 1; 1; 0];
    verty_E1 = [separatrix_pts(:,2); 0; 0; 0];
    patch(vertx_E1, verty_E1, color_E1_basin, 'EdgeColor', 'none');
    plot(separatrix_pts(:,1), separatrix_pts(:,2), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2.0);
else
    warning('No valid separatrix points remained after purification; basins not plotted.');
end

%%% ===================================================================
%%% === 6. PLOT REPRESENTATIVE TRAJECTORIES
%%% ===================================================================
initial_conditions = {
    {[0.025, 0.95], [0.05, 0.9], [0.1, 0.8], [0.1, 0.7], [0.1, 0.6], [0.1, 0.5], [0.9, 0.1], [0.98, 0.05]}; 
    {[0.05, 0.4], [0.02, 0.6], [0.05, 0.3], [0.2, 0.3], [0.4, 0.2], [0.8, 0.05], [0.5, 0.1], [0.9, 0.02]}
};
color_E4_traj = [0.1, 0.4, 0.8];
color_E1_traj = [0.8, 0.3, 0];
for i = 1:length(initial_conditions)
    conditions = initial_conditions{i};
    if i == 1; current_color = color_E4_traj; else; current_color = color_E1_traj; end
    
    for j = 1:length(conditions)
        ic = conditions{j};
        plot(ic(1), ic(2), 'o', 'MarkerFaceColor', current_color, 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
        
        [~, Y] = ode45(ode_function, [0 500], ic, options);
       
        plot(Y(:,1), Y(:,2), '-', 'Color', current_color, 'LineWidth', 2);
        
    end
end

%%% ===================================================================
%%% === 7. PLOT EQUILIBRIUM POINTS & FINAL FORMATTING
%%% ===================================================================
% Stable points (ESS): Filled circles, colored to match their basin
plot(0, 0, 'o', 'MarkerSize', 11, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_E1_traj); 
text(0-0.03, 0-0.03, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 1, 'o', 'MarkerSize', 11, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_E4_traj); 
text(1+0.03, 1+0.03, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
% Unstable points: Solid dark gray circles
color_unstable = [0.3 0.3 0.3]; 
plot(0, 1, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_unstable); 
text(0-0.03, 1+0.03, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_unstable); 
text(1+0.03, 0-0.03, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
% Saddle point: Filled black diamond for clear distinction
plot(p_star, q_star, 's', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'k'); 
text(p_star, q_star-0.02, '$E_5$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
% Final plot adjustments
axis([0 1 0 1]); axis square;
xlabel('Family Strategy, $p$ (Probability of Concentration)','Interpreter','latex', 'FontSize', 14);
ylabel('Manager Strategy, $q$ (Probability of Stewardship)','Interpreter','latex', 'FontSize', 14);
title('Phase Portrait for Bistability under Theorem 1','Interpreter','latex', 'FontSize', 16);
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'FontName', 'Helvetica', 'Layer', 'top');
hold off;
% Export the figure
outputFileName = 'S1_phase_portrait.pdf';
exportgraphics(gca, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);

%%% ===================================================================
%%% === LOCAL FUNCTIONS
%%% ===================================================================
function dydt = replicator_dynamics(~, y, phi, m, n, omega, alpha, beta, V1, V2, V3, V4, f1, f_alpha, f_PM)
    p=y(1); q=y(2); dydt=zeros(2,1);
    U_FC = q*(1-omega)*V1 + (1-q)*(1-omega-phi+m)*V2 + f1;
    U_FD = q*(alpha*(1-omega)+(1-beta)*(1-alpha)*(phi-n))*V3 + (1-q)*alpha*(1-omega-phi+m+n)*V4 + f_alpha;
    U_PS = p*omega*V1 + (1-p)*(omega*V3 + beta*(1-alpha)*(phi-n)*V3) + f_PM;
    U_PA = p*(omega+phi-m)*V2 + (1-p)*(omega+phi-m-n)*V4;
    dydt(1) = p*(1-p)*(U_FC-U_FD);
    dydt(2) = q*(1-q)*(U_PS-U_PA);
end
function J = calculate_jacobian(p, q, phi, m, n, omega, alpha, beta, V1, V2, V3, V4, f1, f_alpha, f_PM)
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