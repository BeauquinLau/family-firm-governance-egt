% --- Script to simulate and plot phase portraits for a dynamic system ---
% This script solves a system of two ODEs for various values of a parameter 'n'
% and plots the resulting trajectories in the p-q phase space.

% --- Parameters & Initialization ---
n_values = [0.123, 0.124, 0.125, 0.126, 0.127, 0.13, 0.135, 0.14, 0.147]; % Define array of n values to be tested
colors = jet(length(n_values)); % Generate a colormap with a distinct color for each n value

% --- Figure and Plotting Loop ---
figure;                 % Create a new figure window
hold on;                % Enable overlaying of subsequent plots

% Loop through each value of the parameter n
for k = 1:length(n_values)
    % Get the current value of n for this iteration
    n = n_values(k);
    
    % Create a function handle for the ODE system with the current n
    ode_function_handle = @(t, y) evolutionaryODE(t, y, n);
    
    % Solve the system of ODEs from the initial point 'B' [0.935, 0.146]
    [~, Y] = ode45(ode_function_handle, [0 1000], [0.935, 0.146]);
    
    % Generate a legend label for the current n value, using TeX for subscript
    legend_label = sprintf('n_{%d} = %.3f', k, n);
    
    % Plot the phase trajectory (p vs. q)
    plot(Y(:,1), Y(:,2), 'o-',...
        'Color', colors(k,:),...
        'LineWidth', 2,...
        'DisplayName', legend_label);
end

% --- Figure Annotations and Formatting ---
% Plot fixed points (E1-E4) and hide them from the legend
plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(0.01, -0.06, 'E_1', 'FontSize', 12, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

plot(0, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(0, 1, 'E_2', 'FontSize', 12, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

plot(1, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(1+0.02, 0, 'E_3', 'FontSize', 12, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

plot(1, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(1, 1, 'E_4', 'FontSize', 12, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Mark the initial point 'B'
text(0.935-0.03, 0.146-0.04, 'B', 'FontSize', 13, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Set uniform figure properties
axis([0 1 0 1]);
axis square;
xlabel('p');
ylabel('q');
grid on;
xticks(0:0.1:1);
yticks(0:0.1:1);
box on;

% Configure and display the legend
lgd = legend('Location', 'northeastoutside');
lgd.Title.String = 'Parameter n';
lgd.FontSize = 10;

% Adjust figure size to accommodate the legend
set(gcf, 'Position', [100 100 1000 600]);
hold off;

% --- Save Figure ---
% Get the current figure handle
fig = gcf;

% Set units to inches for consistent sizing
fig.Units = 'inches';
fig.PaperUnits = 'inches';

% Configure paper size to match the on-screen figure size
pos = fig.Position; % [left, bottom, width, height]
paperWidth = pos(3);
paperHeight = pos(4);
fig.PaperSize = [paperWidth, paperHeight];
fig.PaperPosition = [0, 0, paperWidth, paperHeight];

% Specify output filename
outputFileName = fullfile(outputPath, 'F4_2.eps');

% Save the figure as an EPS file
print(fig, outputFileName, '-depsc');
disp(['Figure saved to: ' outputFileName]);


% --- ODE Function Definition ---
% Defines the system of Ordinary Differential Equations (ODEs)
function dy = evolutionaryODE(~, y, n) % 'n' is passed as a parameter
    % Initialize the derivative vector
    dy = zeros(2,1);

    % Model Parameters (constant for all runs)
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
    
    % The system of ODEs
    % dy(1) corresponds to dp/dt
    dy(1) = y(1)*(1 - y(1))*(y(2)*(1 - omega)*V1 + (1 - y(2)) * (1 - omega - phi + m) * V2 - y(2) * (alpha * (1 - omega) + (1 - alpha) * (phi - n)) * V3 - (1 - y(2)) * alpha * (1 - omega - phi + m + n) * V4 + f_1 - f_2);
    
    % dy(2) corresponds to dq/dt
    dy(2) = y(2)*(1 - y(2))*(y(1)*omega*V1 - y(1)*(omega + phi - m)*V2 + (1 - y(1))*omega*V3 - (1 - y(1))*(omega + phi - m - n)*V4 + f_PM);
end