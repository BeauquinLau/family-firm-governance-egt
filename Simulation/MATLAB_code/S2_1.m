% --- Script to plot a phase portrait with multiple trajectories ---
% This script solves a system of ODEs from various initial conditions
% to visualize the basins of attraction in the p-q phase space.

% --- Initialization ---
figure;                 % Create a new figure window
hold on;                % Enable overlaying of subsequent plots

% Define sets of initial conditions and their corresponding plot styles
% Grouping by style makes the code cleaner and easier to manage.
initial_conditions = {
    % Style 1: Blue squares ('b-s')
    {[0.02 0.1], [0.02 0.2], [0.02 0.3], [0.05 0.15], [0.1 0.4]}, 'b-s';
    
    % Style 2: Red circles ('r-o')
    {[0.02 0.35], [0.02 0.4], [0.02 0.5], [0.2 0.5], [0.7 0.6], [0.75 0.6]}, 'r-o';
    
    % Style 3: Cyan diamonds ('c-d')
    {[0.8 0.6], [0.9 0.6], [0.9 0.63], [0.95 0.65], [0.98 0.65]}, 'c-d';
    
    % Style 4: Green asterisks ('g-*')
    {[0.9 0.65], [0.9 0.8], [0.95 0.7], [0.95 0.9]}, 'g-*'
};

% --- Main Plotting Loop ---
% Loop through each group of initial conditions
for i = 1:size(initial_conditions, 1)
    points = initial_conditions{i, 1}; % Get the list of points for the current style
    style = initial_conditions{i, 2};  % Get the plot style for this group
    
    % Loop through each initial point in the current group
    for j = 1:length(points)
        initial_point = points{j};
        
        % Solve the ODE system from the current initial point
        [~, Y] = ode45(@f2, [0 1000], initial_point);
        
        % Plot the trajectory with the specified style
        % Note: No legend is generated for these individual trajectories
        plot(Y(:,1), Y(:,2), style, 'LineWidth', 2, 'HandleVisibility', 'off');
    end
end

% --- Plot Fixed Points ---
% Plot fixed points (E1-E4) and hide them from the legend
plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(0.01, -0.06, 'E_1', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

plot(0, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(0, 1, 'E_2', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

plot(1, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(1+0.02, 0, 'E_3', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

plot(1, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(1, 1, 'E_4', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% --- Final Figure Formatting ---
axis([0 1 0 1]);
axis('square');
xlabel('p');
ylabel('q');
grid on;
xticks(0:0.1:1);
yticks(0:0.1:1);
box on;
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
outputFileName = fullfile(outputPath, 'F2_1.eps');
% Save the figure as an EPS file
print(fig, outputFileName, '-depsc');
disp(['Figure saved to: ' outputFileName]);

% --- ODE Function Definition ---
% Defines the system of Ordinary Differential Equations (ODEs)
function dy = f2(~, y)
    % Initialize the derivative vector
    dy = zeros(2,1);

    % Model Parameters (These are fixed for this simulation)
    phi = 0.2;
    m = 0.1;
    n = 0.05; % The value of n is fixed in this version of the model
    alpha = 0.5;
    V1 = 3;
    V2 = 5;
    V3 = 7;
    V4 = 9;
    omega = 0.3;
    Fm = 1.09;
    F1 = 2;
    F2 = 1.5;
    
    % The system of ODEs
    % dy(1) corresponds to dp/dt
    dy(1) = y(1)*(1 - y(1))*(y(2)*(1 - omega)*V1 + (1 - y(2)) * (1 - omega - phi + m) * V2 - y(2) * (alpha * (1 - omega) + (1 - alpha) * (phi - n)) * V3 - (1 - y(2)) * alpha * (1 - omega - phi + m + n) * V4 + F1 - F2);
    
    % dy(2) corresponds to dq/dt
    dy(2) = y(2)*(1 - y(2))*(y(1)*omega*V1 - y(1)*(omega + phi - m)*V2 + (1 - y(1))*omega*V3 - (1 - y(1))*(omega + phi - m - n)*V4 + Fm);
end