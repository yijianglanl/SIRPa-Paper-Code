
format long 

warning('off','all')

warning('off','all')

m1 = 12;
m2 = 16;

% Lower bounds for datasets
n1 = 5;
n2 = 5;

load('../Data/SIRP2g.txt');
tdata = SIRP2g(:, 1);
Tdata1 = SIRP2g(:, 2);
Tdata2 = SIRP2g(:, 3);

td = [tdata(n1:m1); tdata(n2:m2)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2)];

lb = [100 100  1e1 10 10 0 0 1e-9 ];
ub = [1e8 1e8 1e7 1e7 1e7 5 10 1e-6];
x = [1e6 1e6  1e5 1e5 1e6 3 5 1e-7];

d=8;

    
   options = optimset('Disp', 'iter', 'TolX', 1e-12, 'TolFun', 1e-12);
x = lsqcurvefit(@(x, td) cost_function(x, tdata, n1, n2, m1, m2, d), x, td, Tdata, lb, ub, options);


    
    cost_value = cost_function(x, tdata, n1, n2,  m1, m2,d);
    Rd = norm(cost_value - Tdata) / norm(Tdata);
    % Calculate the Gaussian norm (L2 norm) for error calculation
    error = cost_value - Tdata;
    gaussian_norm = sqrt(sum(error.^2));
    
    % Normalize the Gaussian norm
    Tdata_norm = norm(Tdata);
    Rd2 = gaussian_norm / Tdata_norm;
    
    %--------------------------------------------------------------------------
    % PARAMETER ESTIMATION - ODE23S AND LSQCURVE FIT
    %--------------------------------------------------------------------------
    disp('____________________________________________________________________________________________________')
    disp('..alpha..Beta...lambda...S2...c1...c2...c3...c4..d1...d2...d3...d4...m0. Rd..')
    disp('____________________________________________________________________________________________________')
    
    % CREATION OF THE VECTOR FOR OUTPUT
    fprintf('%10.10f %10.10f %10.5f %10.10f %10.10f %10.5f %10.10f %10.10f %10.10f %10.5f %10.10f %10.5f %10.5f %10.5f\n', ...
        x(1), x(2), x(3), x(4), x(5),x(6),x(7),x(8), Rd);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
residuals = cost_value - Tdata;
n_total = numel(Tdata);
num_params = numel(x);
AIC = n_total * log(sum(residuals.^2) / n_total) + 2 * num_params;
BIC = n_total * log(sum(residuals.^2) / n_total) + num_params * log(n_total);

% Display the AIC and BIC values
fprintf('AIC: %.2f\n', AIC);
fprintf('BIC: %.2f\n', BIC);
fprintf('Rd: %.2f\n', Rd);

% Subdued color options
color1 = [0.8, 0.2, 0.2];  % Dark Red
color2 = [0.2, 0.6, 0.2];  % Dark Green
             
st1=(1/sqrt(12))*[
14532551.8
26752449.99
37020114.65
65751047.5
40525245.9
79800501.25
120323549.4
207571321.1

];
st2=(1/sqrt(16))*[
15600126.62
9230143.853
36896040.14
41996772.36
30069676.29
31961156.98
26860811.79
12952543.96
0
0
0
0
];

% Subdued color options
color1 = [0.8, 0.2, 0.2];  % Dark Red
color2 = [0.2, 0.6, 0.2];  % Dark Green

% Create figure
figure1 = figure;
set(figure1, 'Position', [100, 100, 800, 600]); % Adjust size for better visibility

% Plot all datasets and model results in one figure
errorbar(tdata(n1:m1), Tdata1(n1:m1)./1e6, st1./1e6, 'o', 'Color', color1, 'MarkerSize', 8, 'MarkerFaceColor', color1, 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 3);
hold on;
plot(tdata(n1:m1), cost_value(1:m1-n1+1)./1e6, '-', 'Color', color1, 'LineWidth', 3);

errorbar(tdata(n2:m2), Tdata2(n2:m2)./1e6, st2./1e6, 'o', 'Color', color2, 'MarkerSize', 8, 'MarkerFaceColor', color2, 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 3);
plot(tdata(n2:m2), cost_value(m1-n1+2:m1-n1+m2-n2+2)./1e6, '-', 'Color', color2, 'LineWidth', 3);

% X and Y labels
xlabel('Time (days)', 'FontSize', 19, 'FontWeight', 'bold');
ylabel('Tumor volume (mm^3)', 'FontSize', 19, 'FontWeight', 'bold');
ylim([0.1e8/1e6, 1.4e9/1e6]);
xlim([7.5,26.4])
grid on;

% Removing legend for magenta lines (xline)
% Keep only data legend
legend({'i.v WT BMDMs+8Gy Data','Fit ', ...
        'i.v SIRP\alpha BMDMs+8Gy Data', 'Fit'}, ...
        'FontSize', 20, 'Location', 'southeast');
    
set(gca, 'FontSize', 19, 'LineWidth', 1.3);




function dydt = g1(t, y, x, d)
    dydt = zeros(2, 1);
    t0 = d / 1.2;
    if (d == 0)
        u = 1;
    else
        u = 2 * (2.558775384581 * t0 + exp(-2.558775384581 * t0) - 1) / ((2.558775384581^2) * t0 * t0);
    end
    
    S1 = exp(-0.005147638197 * d - 0.006577733464 * d * d * u);
    
    dydt(1) = 0.453709396815 * y(1) * (1 - y(1) / 1759812730.672900915146) - 0.000000044133 * y(1) * y(2);
    dydt(2) = 16338.485039659030 - 0.069125930809 * y(2) - 0.000000001991 * y(1) * y(2);
    
    if (t >= 12 && t < 14)
        dydt(1) = dydt(1) * S1;
    elseif (t >= 14)
        dydt(1) = dydt(1) * S1;
    else
        dydt(1) = dydt(1);
    end
end

function dydt = g2(t, y, x, d)
    dydt = zeros(3, 1);
    t0 = d / 1.2;
    if (d == 0)
        u = 1;
    else
        u = 2 * (2.558775384581 * t0 + exp(-2.558775384581 * t0) - 1) / ((2.558775384581^2) * t0 * t0);
    end
    
    S1 = exp(-0.005147638197 * d - 0.006577733464 * d * d * u);
    
    dydt(1) = 0.453709396815 * y(1) * (1 - y(1) / 1.91e9) - 0.000000044133 * y(1) * y(2) - 4.43e-8* y(1) * y(3);
    dydt(2) = 16338.485039659030 - 0.069125930809 * y(2) - 0.000000001991 * y(1) * y(2);
    dydt(3) = 0.000000439623 - 0.394943790649 * y(3) - 0.000000815338 * y(1) * y(3);
    
    if (t >= 12 && t < 14)
        dydt(1) = dydt(1) * S1 - dydt(1) * (1 - S1) *8.2474254899;
    elseif (t >= 14)
        dydt(1) = dydt(1) - dydt(1) * (1 - S1) *19.9;
    else
         dydt(1) = dydt(1);
    end
end
   

function cost = cost_function(x, tdata, n1, n2, m1, m2, d)
    [~, Y1_common] = ode23s(@(t, y) g1(t, y, x, d), tdata(n1:m1), [x(1); x(3)]);
    [~, Y2_common] = ode23s(@(t, y) g2(t, y, x, d), tdata(n2:m2), [x(2); x(4);x(5)]);
   
    % Concatenate the sums into a single vector
    cost = [Y1_common(:, 1) + Y1_common(:, 2);
            Y2_common(:, 1) + Y2_common(:, 2)+ Y2_common(:, 3) ;
        ];
    

    % Ensure that the size of the cost matches the size of Tdata
    assert(size(cost, 1) == size(tdata(n1:m1), 1) + size(tdata(n2:m2), 1) , 'Inconsistent dimensions');
end