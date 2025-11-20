
format long 

warning('off','all')

warning('off','all')
load('../Data/SIRP2e.txt');
% Upper bounds for datasets
m1 = 12;
m2 = 14;
m3 = 14;
m4 = 14;
% Lower bounds for datasets
n1 = 5;
n2 = 5;
n3 = 5;
n4 = 5;

tdata = SIRP2e(:, 1);
Tdata1 = SIRP2e(:, 2);
Tdata2 = SIRP2e(:, 3);
Tdata3 = SIRP2e(:, 4);
Tdata4 = SIRP2e(:, 5);


lb = [100 100 100  100 1e1 1e1 1e1 0];
ub = [Tdata1(n1+1)  Tdata2(n2+1) Tdata3(n3+1) Tdata4(n4+1) 1e7 1e7 1e7 10];
x = [Tdata1(n1)  Tdata2(n2) Tdata3(n3) Tdata4(n4) 1e5 1e5 1e5 5];

d1=0;
d2=8;

td = [tdata(n1:m1); tdata(n2:m2); tdata(n3:m3); tdata(n4:m4)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2); Tdata3(n3:m3); Tdata4(n4:m4)];

    
   options = optimset('Disp', 'iter', 'TolX', 1e-12, 'TolFun', 1e-12);
x = lsqcurvefit(@(x, td) cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2]), x, td, Tdata, lb, ub, options);


    
    cost_value = cost_function(x, tdata, n1, n2, n3,n4, m1, m2, m3,m4,[d1, d2]);
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
        x(1), x(2), x(3), x(4), x(5), x(6), x(7), Rd);
    
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

% color options
color1 = [0.8, 0.2, 0.2];  
color2 = [0.2, 0.6, 0.2];  
color3 = [0.2, 0.2, 0.6];  
color4 = [0.2, 0.2, 0.2];                   
st1=(1/sqrt(12))*[
19394831.94
19413937.97
29046352.91
94614916.74
132711333.4
188509516.4
227539502.1
122600000
];
st2=(1/sqrt(14))*[
8980775.077
11537475.17
23562823.83
22775999.93
43456622.33
50618510.25
33963491.66
10339472.49
0
0
];
st3=(1/sqrt(14))*[
9434805.539
24362037.86
31648099.82
38680141.33
98530590.6
153496561.3
228939191.3
314508290.3
139244260.5
160186231.9
];
st4=(1/sqrt(14))*[
9825904.287
19471371.83
27783359.44
70585062.95
146576532.9
172911159.7
212891238.7
188493908.6
201953713.5
0
];
 figure1 = figure;
set(figure1, 'Position', [100, 100, 1200, 600]); 

% Plot all datasets and model results 
plot(tdata(n1:m1), cost_value(1:m1-n1+1)./1e6, '-', 'Color', color1, 'LineWidth', 3);hold on;
errorbar(tdata(n1:m1), Tdata1(n1:m1)./1e6, st1./1e6, 'o', 'Color', color1, 'MarkerSize', 8, 'MarkerFaceColor', color1, 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 3);
hold on;

plot(tdata(n2:m2), cost_value(m1-n1+2:m1-n1+m2-n2+2)./1e6, '-', 'Color', color2, 'LineWidth', 3);
errorbar(tdata(n2:m2), Tdata2(n2:m2)./1e6, st2./1e6, 'o', 'Color', color2, 'MarkerSize', 8, 'MarkerFaceColor', color2, 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 3);

plot(tdata(n3:m3), cost_value(m1-n1+m2-n2+3:m1-n1+m2-n2+m3-n3+3)./1e6, '-', 'Color', color3, 'LineWidth', 3);
errorbar(tdata(n3:m3), Tdata3(n3:m3)./1e6, st3./1e6, 'o', 'Color', color3, 'MarkerSize', 8, 'MarkerFaceColor', color3, 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 3);

plot(tdata(n4:m4), cost_value(m1-n1+m2-n2+m3-n3+4:m1-n1+m2-n2+m3-n3+m4-n4+4)./1e6, '-', 'Color', color4, 'LineWidth', 3);
errorbar(tdata(n4:m4), Tdata4(n4:m4)./1e6, st4./1e6, 'o', 'Color', color4, 'MarkerSize', 8, 'MarkerFaceColor', color4, 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 3);

%xline(12,'--', 'm', 'LineWidth', 2, 'Label', 'RT Treatment','fontsize',19, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'middle');
xline(12, '--m', 'LineWidth', 2, 'Label', '8Gy Treatment','fontsize',19, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'middle');
xlim([7.6, 26.5]);

xlabel('Time (days)', 'FontSize', 19, 'FontWeight', 'bold');
ylabel('Tumor volume (mm^3)', 'FontSize', 19, 'FontWeight', 'bold');
ylim([0.1e8/1e6, 2e9/1e6]);
grid on;
legend({'Prediction','SIRP\alpha+0R Data', 'Prediction', 'Cl2MDP+8Gy Data', 'Prediction', ...
        'αCSF1R+8GY Data', 'Prediction', 'SIRP\alpha+8Gy Data'}, ...
        'FontSize', 20, 'Location', 'northwest');
set(gca, 'FontSize', 19, 'LineWidth', 1.3);



    function dydt = g(t, y, x, d)
    dydt = zeros(2,1);
    t0 = d / 1.2;
    if d == 0
        u = 1;
    else
        u = 2 * (2.558775384581 * t0 + exp(-2.558775384581 * t0) - 1) / ((2.558775384581^2) * t0 * t0);
    end
    
    S1 = exp(-0.005147638197 * d - 0.006577733464 * d * d * u);
    dydt(1) = 0.453709396815 * y(1) * (1 - y(1) / 1759812730.672900915146) - 0.000000044133 * y(1) * y(2);
    dydt(2) = 16338.485039659030 - 0.069125930809 * y(2) - 0.000000001991 * y(1) * y(2);
    
    if t >= 12
        dydt(1) = dydt(1) * S1;
    end
end

function dydt = g1(t, y, x, d)
    dydt = zeros(3,1);
    t0 = d / 1.2;
    if d == 0
        u = 1;
    else
        u = 2 * (2.558775384581 * t0 + exp(-2.558775384581 * t0) - 1) / ((2.558775384581^2) * t0 * t0);
    end
    
    S1 = exp(-0.005147638197 * d - 0.006577733464 * d * d * u);
    dydt(1) = 0.453709396815 * y(1) * (1 - y(1) / 1759812730.672900915146) - 0.000000044133 * y(1) * y(2) - 4.43e-8 * y(1) * y(3);
    dydt(2) = 16338.485039659030 - 0.069125930809 * y(2) - 0.000000001991 * y(1) * y(2);
    dydt(3) = 0.000000439623 - 0.394943790649 * y(3) - 0.000000815338 * y(1) * y(3);
    A=15.3;
    if t >= 12
        dydt(1) = dydt(1) * S1 - dydt(1) * (1 - S1) * A;
    end
end
   

function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, d)
    [~, Y1_common] = ode23s(@(t, y) g1(t, y, x, d(1)), tdata(n1:m1), [x(1); x(5);x(6)]);
    [~, Y2_common] = ode23s(@(t, y) g1(t, y, x, d(2)), tdata(n2:m2), [x(2); x(5);x(6)]);
    [~, Y3_common] = ode23s(@(t, y) g(t, y, x, d(2)), tdata(n3:m3), [x(3);x(7)]);
    [~, Y4_common] = ode23s(@(t, y) g(t, y, x, d(2)), tdata(n4:m4), [x(4);x(7)]);

    % Concatenate the sums into a single vector
    cost = [Y1_common(:, 1) + Y1_common(:, 2)+ Y1_common(:, 3) ;
            Y2_common(:, 1) + Y2_common(:, 2)+ Y2_common(:, 3) ;
            Y3_common(:, 1) + Y3_common(:, 2) ;
            Y4_common(:, 1) + Y4_common(:, 2) ];
    

    % size of the cost matches the size of Tdata
    assert(size(cost, 1) == size(tdata(n1:m1), 1) + size(tdata(n2:m2), 1) + size(tdata(n3:m3), 1) + size(tdata(n4:m4), 1), 'Inconsistent dimensions');
end
