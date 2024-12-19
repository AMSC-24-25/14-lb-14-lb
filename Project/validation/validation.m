% Specify filename
fileName = 'simulation_data.csv';

% Read from CSV files
data = readmatrix(fileName);
u_paper = readmatrix('u_velocity_center_vertical_paper.csv');
v_paper = readmatrix('v_velocity_center_horizontal_paper.csv');

% Extract coordinates and mesaurments
X = data(data(:,3) == 64, 2);
v = data(data(:,3) == 64, 6) ./ (100/(6*129)); 

Y = data(data(:,2) == 64, 3);
u = data(data(:,2) == 64, 5) ./ (100/(6*129));

% Plot v
% Create plot
figure;
plot(v_paper(:,1), v_paper(:,3), 'x', 'DisplayName', 'v paper');
hold on;

plot(X, v, '-', 'DisplayName', 'v LBM');

xlabel('Horizontal Position');
ylabel('v magnitude normalized to lid velocity');
title('Comparison between Ghia 1982 and LBM implementation');
legend;
grid on;

hold off;

% Plot u
% Create plot
figure;
plot(u_paper(:,2), u_paper(:,3), 'x', 'DisplayName', 'u paper');
hold on;

plot(Y, u, '-', 'DisplayName', 'u LBM');

xlabel('Vertical Position');
ylabel('u magnitude normalized to lid velocity');
title('Comparison between Ghia 1982 and LBM implementation');
legend;
grid on;

hold off;
