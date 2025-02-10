clear; close all;

% Load data
load('dataset_A_B.mat');  

% Number of matrices (num_matrices = 50 in manuscript)
num_matrices = 5;

% Initialize cell arrays (3 methods, 50 matrices each)
time_tracking = cell(3, num_matrices);
objerror_tracking = cell(3, num_matrices);
solerror_tracking = cell(3, num_matrices);

% Loop through matrices
for k = 1:num_matrices
    % Extract current A, B, and P
    A_k = eval(['A' num2str(k)]);
    B_k = eval(['B' num2str(k)]);
    P_k = eval(['P' num2str(k)]);
    
    % Method 1: project_reweighted
    tic;
    [X2, f_best2, outref] = linear_reweighted_solver(P_k, A_k, B_k, 1, 10, 1e-5);
    time_tracking{1, k} = outref.time_tracking;
    objerror_tracking{1, k} = outref.objerror_tracking;
    solerror_tracking{1, k} = outref.solerror_tracking;

    % Method 2: lp_regularizationt (p = 0.75)
    tic;
    [X3, f_best3, outlp75] = lp_norm_solver(P_k, A_k, B_k, 0.75, 1, 10, 1e-5);
    time_tracking{2, k} = outlp75.time_tracking;  
    objerror_tracking{2, k} = outlp75.objerror_tracking;
    solerror_tracking{2, k} = outlp75.solerror_tracking;

    % Method 3: lp_regularizationt (p = 0.5)
    tic;
    [X4, f_best4, outlp5] = lp_norm_solver(P_k, A_k, B_k, 0.5, 1, 10, 1e-5);
    time_tracking{3, k} = outlp5.time_tracking;
    objerror_tracking{3, k} = outlp5.objerror_tracking;
    solerror_tracking{3, k} = outlp5.solerror_tracking;
end

% Save results
save('tracking_data02.mat', 'time_tracking', 'objerror_tracking', 'solerror_tracking');

%% Load Data for Plotting
load('tracking_data02.mat');
num_matrices = 5;
solver_names = {'linear reweighted', '$L_{p=0.75}$ norm', '$L_{p=0.5}$ norm'};
colors = [1, 0, 0;        % reweighted (Red)
          0.1, 0.1, 0.8;  % lp_norm_0.75 (Blue)
          0.8, 0.8, 0.1]; % lp_norm_0.5 (Yellow)];

target_times = 0:0.1:4; % Define interpolation times

%% **Plot 1: Objective Error vs. Time (Log Scale)**
figure;
hold on;
plot_handles = gobjects(1, 3);

for i = 1:3
    new_error_vals = zeros(num_matrices, length(target_times));

    for k = 1:num_matrices
        time_vals = time_tracking{i, k};
        error_vals = objerror_tracking{i, k};

        % Interpolate error values at target times
        for q = 1:length(target_times)
            target_time = target_times(q);
            [~, closest_index] = min(abs(time_vals - target_time));
            new_error_vals(k, q) = error_vals(closest_index);
        end
    end

    % Compute mean and standard deviation
    averagevalue = mean(new_error_vals, 1);
    sdvalue = std(new_error_vals, 0, 1);

    % Plot shaded area for standard deviation
    fill([target_times, fliplr(target_times)], ...
    [averagevalue + sdvalue, fliplr(max(averagevalue - sdvalue, 1e-2*ones(size(averagevalue))))], ...
    colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');  

    % Plot mean line
    plot_handles(i) = plot(target_times, averagevalue, 'o-', ...
                           'LineWidth', 2, 'MarkerSize', 4, ...
                           'Color', colors(i, :));
end

% Configure plot
xlim([0, 3]);
set(gca, 'YScale', 'log'); % Set log scale for y-axis
xlabel('Time (seconds)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Objective Error', 'Interpreter', 'latex', 'FontSize', 20);
title('Objective Error vs. Time', 'Interpreter', 'latex', 'FontSize', 20);
legend(plot_handles, solver_names, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);
grid on;
hold off;

%% **Plot 2: Residual vs. Time (Log Scale)**
figure;
hold on;
plot_handles2 = gobjects(1, 3);

for i = 1:3
    new_error_vals = zeros(num_matrices, length(target_times));

    for k = 1:num_matrices
        time_vals = time_tracking{i, k};
        error_vals = solerror_tracking{i, k}; % Residual error

        % Interpolate error values at target times
        for q = 1:length(target_times)
            target_time = target_times(q);
            [~, closest_index] = min(abs(time_vals - target_time));
            new_error_vals(k, q) = error_vals(closest_index);
        end
    end

    % Compute mean and standard deviation
    averagevalue = mean(new_error_vals, 1);
    sdvalue = std(new_error_vals, 0, 1);

    % Plot shaded area for standard deviation
    fill([target_times, fliplr(target_times)], ...
         [averagevalue + sdvalue, fliplr(max(averagevalue - sdvalue, 1e-1*ones(size(averagevalue))))], ...
         colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 

    % Plot mean line
    plot_handles2(i) = plot(target_times, averagevalue, 'o-', ...
                            'LineWidth', 2, 'MarkerSize', 4, ...
                            'Color', colors(i, :));
end

% Configure plot
xlim([0, 3]);
set(gca, 'YScale', 'log'); % Set log scale for y-axis
xlabel('Time (seconds)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Residual', 'Interpreter', 'latex', 'FontSize', 20);
title('Residual vs. Time', 'Interpreter', 'latex', 'FontSize', 20);
legend(plot_handles2, solver_names, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);
grid on;
hold off;


