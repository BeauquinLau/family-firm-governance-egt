% --- Script to simulate and plot the time evolution of a dynamic system ---
% This script solves a system of two ODEs for various values of a parameter 'n'
% and plots the time evolution of the first state variable, p.

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
    
    % Solve the system of ODEs
    % The initial condition is [p(0), q(0)] = [0.935, 0.146]
    [T, Y] = ode45(@(t, y) dynamic_system(t, y, n), [0 1000], [0.935, 0.146]);
    
    % Generate a legend label for the current n value, using TeX for subscript
    legend_label = sprintf('n_{%d} = %.3f', k, n);
    
    % Plot the time evolution of p (the first state variable, Y(:,1))
    plot(T, Y(:,1), '-',...
        'Color', colors(k,:),...
        'LineWidth', 2,...
        'DisplayName', legend_label);
end

% --- Figure Formatting ---
xlim([0 1000]);      % Set the x-axis (time) limits
ylim([0 1]);         % Set the y-axis (probability p) limits
xlabel('Time (t)');
ylabel('p');
grid on;

% Configure and display the legend
lgd = legend('Location', 'northeastoutside');
lgd.Title.String = 'Parameter n';
lgd.FontSize = 10;

% Adjust figure size to accommodate the legend
set(gcf, 'Position', [100 100 1000 600]);
box on;
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
outputFileName = fullfile(outputPath, 'F6_1.eps');

% Save the figure as an EPS file
print(fig, outputFileName, '-depsc');
disp(['Figure saved to: ' outputFileName]);


% --- ODE Function Definition ---
% Defines the system of Ordinary Differential Equations (ODEs)
function dy = dynamic_system(~, y, n) % 'n' is passed as a parameter
    % Initialize the derivative vector
    dy = zeros(2, 1);
    
    % Model Parameters (constant for all runs)
    phi = 0.25;
    m = 0.05;
    alpha = 0.5;
    V1 = 1;
    V2 = 2;
    V3 = 3;
    V4 = 5;
    omega = 0.2;
    Fm = 0.63;
    F1 = 1.6;
    F2 = 1;
    
    % The system of ODEs
    % dy(1) corresponds to dp/dt
    dy(1) = y(1) * (1 - y(1)) * (y(2) * (1 - omega) * V1 + ...
        (1 - y(2)) * (1 - omega - phi + m) * V2 - ...
        y(2) * (alpha * (1 - omega) + (1 - alpha) * (phi - n)) * V3 - ...
        (1 - y(2)) * alpha * (1 - omega - phi + m + n) * V4 + F1 - F2);
    
    % dy(2) corresponds to dq/dt
    dy(2) = y(2) * (1 - y(2)) * (y(1) * omega * V1 - ...
        y(1) * (omega + phi - m) * V2 + ...
        (1 - y(1)) * omega * V3 - ...
        (1 - y(1)) * (omega + phi - m - n) * V4 + Fm);
end