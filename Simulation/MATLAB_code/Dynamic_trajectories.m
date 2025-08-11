% --- Script to simulate and plot stable/unstable manifolds in a phase portrait ---
% This script calculates saddle points for various values of a parameter 'n',
% linearizes the system at these points, and plots the resulting manifolds.

% --- Model Parameters ---
phi = 0.25;
m = 0.05;
alpha = 0.5;
V1 = 1;
V2 = 2;
V3 = 3;
V4 = 5;
omega = 0.2;
f_PM = 0.63;
f_1 = 1.6;
f_2 = 1;

% --- Simulation Setup ---
n_values = [0.123,0.124, 0.125, 0.126, 0.127, 0.13, 0.135, 0.14, 0.147]; % Define array of n values to be tested
p_star = zeros(size(n_values)); % Pre-allocate array for p* results
q_star = zeros(size(n_values)); % Pre-allocate array for q* results

% --- Calculate Equilibrium Points (Saddle Points) ---
% Loop through each n to find the corresponding saddle point (p*, q*)
for i = 1:length(n_values)
    n = n_values(i);
    
    % Calculate p*
    numerator_p = (omega + phi - m - n) * V4 - omega * V3 - f_PM;
    denominator_p = omega * V1 - (omega + phi - m) * V2 - omega * V3 + (omega + phi - m - n) * V4;
    p_star(i) = numerator_p / denominator_p;
    
    % Calculate q*
    numerator_q = alpha * (1 - omega - phi + m + n) * V4 - (1 - omega - phi + m) * V2 - f_1 + f_2;
    denominator_q = (1 - omega) * V1 - (1 - omega - phi + m) * V2 - (alpha * (1 - omega) + (1 - alpha) * (phi - n)) * V3 + alpha * (1 - omega - phi + m + n) * V4;
    q_star(i) = numerator_q / denominator_q;
end

% --- Integration & Figure Setup ---
eps = 1e-5; % Perturbation size for manifold calculation
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9); % Set ODE solver options

figure;
hold on;
colors = jet(length(n_values)); % Generate a colormap for plotting
legend_handles = gobjects(length(n_values),1); % Pre-allocate array for legend handles

% --- Main Loop: Plot Manifolds for Each Saddle Point ---
for i = 1:length(n_values)
    n = n_values(i);
    y1_val = p_star(i);
    y2_val = q_star(i);
    
    % --- Jacobian Calculation (Numerical) ---
    y1 = y1_val;
    y2 = y2_val;
    G_y2 = y2*(1 - omega)*V1 + (1 - y2)*(1 - omega - phi + m)*V2 - y2*(alpha*(1 - omega) + (1 - alpha)*(phi - n))*V3 - (1 - y2)*alpha*(1 - omega - phi + m + n)*V4 + f_1 - f_2;
    H_y1 = y1*omega*V1 - y1*(omega + phi - m)*V2 + (1 - y1)*omega*V3 - (1 - y1)*(omega + phi - m - n)*V4 + f_PM;
    dG_dy2 = (1 - omega)*V1 - (1 - omega - phi + m)*V2 - (alpha*(1 - omega) + (1 - alpha)*(phi - n))*V3 + alpha*(1 - omega - phi + m + n)*V4;
    dH_dy1 = omega*V1 - (omega + phi - m)*V2 - omega*V3 + (omega + phi - m - n)*V4;
    j11 = (1 - 2*y1) * G_y2;
    j12 = y1 * (1 - y1) * dG_dy2;
    j21 = y2 * (1 - y2) * dH_dy1;
    j22 = (1 - 2*y2) * H_y1;
    J_E = [j11, j12; j21, j22]; % Assemble the Jacobian matrix
    
    % --- Eigenvalue/Eigenvector Analysis ---
    [V, D] = eig(J_E);
    lambda = diag(D);
    [~, idx_pos] = max(real(lambda)); % Index of unstable eigenvector
    [~, idx_neg] = min(real(lambda)); % Index of stable eigenvector
    v_unstable = V(:, idx_pos);
    v_stable = V(:, idx_neg);
    
    % --- Manifold Integration ---
    % Integrate along the unstable manifold (forward in time)
    init_unstable_pos = [y1_val; y2_val] + eps * v_unstable;
    init_unstable_neg = [y1_val; y2_val] - eps * v_unstable;
    [~, Y_uni_pos] = ode45(@(t, y) f1(t, y, n), [0 1000], init_unstable_pos, options);
    [~, Y_uni_neg] = ode45(@(t, y) f1(t, y, n), [0 1000], init_unstable_neg, options);
    
    % Integrate along the stable manifold (backward in time)
    init_stable_pos = [y1_val; y2_val] + eps * v_stable;
    init_stable_neg = [y1_val; y2_val] - eps * v_stable;
    [~, Y_sta_pos] = ode45(@(t, y) f1(t, y, n), [0 -1000], init_stable_pos, options);
    [~, Y_sta_neg] = ode45(@(t, y) f1(t, y, n), [0 -1000], init_stable_neg, options);
    
    % --- Plotting ---
    % Plot the primary trajectory (one branch of the unstable manifold) and save handle for legend
    h = plot(Y_uni_pos(:,1), Y_uni_pos(:,2), '--', 'Color', colors(i,:),...
        'LineWidth', 2, 'DisplayName', sprintf('n_{%d}=%.3f', i, n_values(i))); 
    legend_handles(i) = h;
    
    % Plot the other three manifold branches without legend entries
    plot(Y_uni_neg(:,1), Y_uni_neg(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off'); 
    plot(Y_sta_pos(:,1), Y_sta_pos(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    plot(Y_sta_neg(:,1), Y_sta_neg(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
 
    % Mark the saddle point, hide from legend
    plot(y1_val, y2_val, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','k', 'HandleVisibility','off'); 
end

% Plot other fixed points and markers, hidden from legend
plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(0, 0, 'E_1', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(0, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(0, 1, 'E_2', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(1, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(1+0.02, 0, 'E_3', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(1, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(1, 1, 'E_4', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(0.4, 0.3, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
text(0.402, 0.302, 'A', 'FontSize', 13, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(0.935, 0.146, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
text(0.937, 0.148, 'B', 'FontSize', 13, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(0.8, 0.85, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
text(0.802, 0.852, 'C', 'FontSize', 13, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% --- Final Figure Formatting ---
axis([0 1 0 1]);
axis('square');
xlabel('p');
ylabel('q');
grid on;
box on;

% Configure and display the legend with TeX interpreter
lgd = legend(legend_handles, 'Interpreter','tex'); 
lgd.Title.String = 'Parameter n';  
set(lgd, 'Location','northeastoutside');

% Adjust figure size
set(gcf, 'Position', [100 100 1000 600]);
hold off;

% --- Save Figure ---
% Get the current figure handle
fig = gcf;
% Set units to inches for consistent sizing
fig.Units = 'inches';
fig.PaperUnits = 'inches';
% Configure paper size to match the on-screen figure size
pos = fig.Position;
paperWidth = pos(3);
paperHeight = pos(4);
fig.PaperSize = [paperWidth, paperHeight];
fig.PaperPosition = [0, 0, paperWidth, paperHeight];
% Specify output filename
outputFileName = fullfile(outputPath, 'F3.eps');
% Save the figure as an EPS file
print(fig, outputFileName, '-depsc');
disp(['Figure saved to: ' outputFileName]);

% --- ODE Function Definition ---
% Defines the system of Ordinary Differential Equations (ODEs)
function dydt = f1(~, y, n)
    % Model Parameters (constant for all runs)
    phi = 0.25; m = 0.05; alpha = 0.5;
    V1 = 1; V2 = 2; V3 = 3; V4 = 5; omega = 0.2;
    f_PM = 0.63; f_1 = 1.6; f_2 = 1;
    % State variables
    y1 = y(1);
    y2 = y(2);
    % The system of ODEs
    dy1dt = y1*(1 - y1)*(y2*(1 - omega)*V1 + (1 - y2)*(1 - omega - phi + m)*V2 ...
        - y2*(alpha*(1 - omega) + (1 - alpha)*(phi - n))*V3 ...
        - (1 - y2)*alpha*(1 - omega - phi + m + n)*V4 + f_1 - f_2);
    dy2dt = y2*(1 - y2)*(y1*omega*V1 - y1*(omega + phi - m)*V2 ...
        + (1 - y1)*omega*V3 - (1 - y1)*(omega + phi - m - n)*V4 + f_PM);
    dydt = [dy1dt; dy2dt];
end