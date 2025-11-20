function model_uncertainty_analysis
    clear;
    close all;
    clc;
    format long;
    warning('off', 'all');
clc;
clf;
clear all;
clearvars -global


mf = mfilename('fullpath');
if isempty(mf)
    
    project_root = pwd;
else
    project_root = fileparts(mf);
end


paths.models    = fullfile(project_root, 'Models');
addpath(genpath(paths.models));

addpath('../Models');
 
m3 = 11;

n3 = 4;


load('../Data/MCS_cells.txt');

tdata = MCS_cells(:,1);

Tdata3 = MCS_cells(:,4);

    % Define the subset of data to use
    td = tdata(n3:m3);
    Tdata = Tdata3(n3:m3);
    % Set parameter bounds and initial guess
    lb = [1e-7, 1e-8, 2e-3, 1e-8,  1e3, 1e5, 1e5];
    ub = [1e-5, 1e-6, 8e-1, 1e-6,  1e8, 1e7, 1e7];
    x0 = [1e-6, 1e-7, 2e-1 , 1e-7, 1e6, 1e6,1e6];

    % Fit model using nonlinear least squares
    options = optimset('Display', 'iter', 'TolX', 1e-15, 'TolFun', 1e-15);
    x = lsqcurvefit(@(x,td) cost_function(x,td,Tdata), x0, td, Tdata, lb, ub, options);

    % Generate model predictions for original fit
    [~, modelFit] = ode23s(@(t, y) MC38_SIRPA_KO(t, y, x), tdata, [x(5); x(6);x(7)]);

    % Define the number of bootstrap iterations and preallocate
    num_iterations =500;
    bootstrap_params = zeros(length(x), num_iterations);

    % Bootstrap process
    for iter = 1:num_iterations
        % Generate Poisson counts and create bootstrap sample
        counts = poissrnd(1, length(Tdata), 1);
        bootstrap_td = repelem(td, counts);
        bootstrap_Tdata = repelem(Tdata, counts);

        if numel(unique(bootstrap_td)) < 2
            continue;
        end

        try
            x_bootstrap = lsqcurvefit(@(x, td) cost_function(x, td, bootstrap_Tdata), x0, bootstrap_td, bootstrap_Tdata, lb, ub, options);
            bootstrap_params(:, iter) = x_bootstrap;
        catch
            continue;
        end
    end

     % Compute confidence intervals
    confidence_intervals = prctile(bootstrap_params', [5, 95]);
    % Output parameter estimates and confidence intervals 
fprintf('Parameter Estimates and Confidence Intervals:\n');
 for i = 1:length(x) 
     fprintf('Parameter %.d: Mean = %.12f, CI = [%.12f, %.12f]\n', i, mean(bootstrap_params(i, :)), confidence_intervals(1, i), confidence_intervals(2, i));
 end
    ci_params_lower = prctile(bootstrap_params', 5);
    ci_params_upper = prctile(bootstrap_params', 95);

    % Model predictions for confidence intervals
    [~, ci_modelFit_lower] = ode23s(@(t, y) MC38_SIRPA_KO(t, y, ci_params_lower), tdata, [ci_params_lower(5); ci_params_lower(6); ci_params_lower(7)]);
    [~, ci_modelFit_upper] = ode23s(@(t, y) MC38_SIRPA_KO(t, y, ci_params_upper), tdata, [ci_params_upper(5); ci_params_upper(6); ci_params_lower(7)]);
    
    % Plotting
    % Ensure that all segments are column vectors
tdata_segment = tdata(n3:m3);
lower_CI_segment = ci_modelFit_lower(1:m3-n3+1);
upper_CI_segment = ci_modelFit_upper(1:m3-n3+1);

% Make sure all are column vectors for consistency
tdata_segment = tdata_segment(:);
lower_CI_segment = lower_CI_segment(:);
upper_CI_segment = upper_CI_segment(:);

% Validate the sizes match
assert(length(tdata_segment) == length(lower_CI_segment) && length(tdata_segment) == length(upper_CI_segment), 'Vector lengths do not match.');

% Plotting setup
figure;
hold on;

% Original Data Points
plot(tdata_segment, Tdata3(n3:m3), 'k*','Markersize',10, 'DisplayName', 'Original Data');

% Model Fit Line
plot(tdata_segment, modelFit(1:length(tdata_segment)), 'b-', 'LineWidth', 2, 'DisplayName', 'Model Fit');

% Confidence Interval Lines
plot(tdata_segment, lower_CI_segment, ':r', 'LineWidth', 1.5, 'DisplayName', 'Lower 5% CI');
plot(tdata_segment, upper_CI_segment, ':r', 'LineWidth', 1.5, 'DisplayName', 'Upper 95% CI');

% Filling between Confidence Intervals
fill([tdata_segment; flipud(tdata_segment)], [lower_CI_segment; flipud(upper_CI_segment)],[0.3, 0.3, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Labels, Title, and Legend
xlabel('Time');
ylabel('Tumor Cells Population');
%title('Model Fit with 95% CI');
legend('Data','Model fit','Location', 'northwest');

% Final adjustments
set(gcf, 'Color', 'w'); 
hold off;
 
end

function cost = cost_function(x, td, Tdata)
    y0 = [x(5); x(6);x(7)];
    [~, Y] = ode23s(@(t, y) MC38_SIRPA_KO(t, y, x), td, y0);
    cost = Y(:, 1)+Y(:,2)+Y(:,3) ;  % Compute residuals for lsqcurvefit
end