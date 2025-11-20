format long 
warning('off','all')

m1 = 16; m2 = 14; m3 = 11; m4 = 8;
n1 = 8;  n2 = 7;  n3 = 4;  n4 = 2;

% ---------------- Data ----------------
load('../Data/WT.txt');

tdata  = WT(:,1);
Tdata1 = WT(:,2);
Tdata2 = WT(:,3);
Tdata3 = WT(:,4);
Tdata4 = WT(:,5);

td    = [tdata(n1:m1); tdata(n2:m2); tdata(n3:m3); tdata(n4:m4)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2); Tdata3(n3:m3); Tdata4(n4:m4)];

% ---------------- Bounds, initial guess ----------------
lb = [0.2 5e8   1e-8   1.3e3  2e-3   3.34e-11  1    1    1    1    1   ];
ub = [0.6 5e9   1e-6   1.3e5  2e-1   3.34e-9   Tdata1(n1+1) Tdata2(n2+1) Tdata3(n3+1) Tdata4(n4+1) 1e7 ];
x0 = [0.45 2e9  1e-7   1.3e4  2e-2   3.34e-10  Tdata1(n1)   Tdata2(n2)   Tdata3(n3)   Tdata4(n4)   1e6];

% ---------------- Parameter estimation ----------------
x = lsqcurvefit(@(x,td) cost_function(x,tdata,n1,n2,n3,n4,m1,m2,m3,m4), ...
                x0, td, Tdata, lb, ub, ...
                optimset('Disp','iter','TolX',1e-15,'TolFun',1e-15));

cost_value = cost_function(x,tdata,n1,n2,n3,n4,m1,m2,m3,m4);

% ---------------- Sensitivity analysis ----------------
sensitivity_on = true;
if sensitivity_on
    
    % Number of parameters to investigate (first 6 here)
    n_par = 6;
    sens_mat = zeros(length(td), n_par); 
    
    % Baseline model output
    y0 = cost_value;

    % Perturbation factor
    perturbation_factor = 0.25;

    for i = 1:n_par
        x_perturbed = x;
        delta = x_perturbed(i) * perturbation_factor;
        x_perturbed(i) = x_perturbed(i) + delta;

        % Evaluate the cost function for the perturbed parameter
        y1 = cost_function(x_perturbed, tdata, n1, n2, n3, n4, m1, m2, m3, m4);

        % Relative sensitivity at each time point
        sens_mat(:, i) = (y1 - y0) ./ max(1e-12, abs(y0));
    end

    % Average absolute sensitivity over all time points
    average_sensitivity = mean(abs(sens_mat), 1);

    % ---------------- Plot: Average sensitivity bar chart ----------------
    figure;
    numColors = n_par;
    colors = jet(numColors);

    hold on;
    for i = 1:numColors
        bar(i, average_sensitivity(i), 'FaceColor', colors(i,:), 'BarWidth',0.5);
    end
    hold off;

    xticks(1:numColors);
    xticklabels({'c_1', 'c_{max}', '\phi_e', '\gamma_e', '\delta_e', '\eta_e'});
    xtickangle(45);

    xlabel('Parameters', 'FontSize', 23);
    ylabel('Sensitivity Index', 'FontSize', 23);
    set(gca, 'FontSize', 18, 'Box','on','BoxStyle','back');
    ylim([0, max(average_sensitivity)*1.2]);
    
end

function dydt = tumor_growth(~,y,x)

dydt = zeros(2,1);
dydt(1) = x(1)*y(1)*(1 - y(1)/x(2)) - x(3)*y(1)*y(2);
dydt(2) = x(4) - x(5)*y(2) - x(6)*y(1)*y(2);

end

% ---------------- Local cost function ----------------
function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4)
    [~,Y1] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n1:m1), [x(7);  x(11)]);
    [~,Y2] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n2:m2), [x(8);  x(11)]);
    [~,Y3] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n3:m3), [x(9);  x(11)]);
    [~,Y4] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n4:m4), [x(10); x(11)]);

    Y1_sum = Y1(:,1) + Y1(:,2);
    Y2_sum = Y2(:,1) + Y2(:,2);
    Y3_sum = Y3(:,1) + Y3(:,2);
    Y4_sum = Y4(:,1) + Y4(:,2);

    cost = [Y1_sum; Y2_sum; Y3_sum; Y4_sum];
end
