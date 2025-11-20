% ====== Calibrate step 3  ======
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
format long
warning('off','all')

% Upper bounds for datasets
m1 = 12; m2 = 12; m3 = 13; m4 = 13;
% Lower bounds for datasets
n1 = 4; n2 = 5;  n3 = 5;  n4 = 5;

load('../Data/WTl100.txt');

tdata  = WTl100(:, 1);
Tdata1 = WTl100(:, 2);
Tdata2 = WTl100(:, 3);
Tdata3 = WTl100(:, 4);
Tdata4 = WTl100(:, 5);

td    = [tdata(n1:m1); tdata(n2:m2); tdata(n3:m3); tdata(n4:m4)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2); Tdata3(n3:m3); Tdata4(n4:m4)];

lb = [5.15e-3  0.0001  0.5   100   100    100         100          1e5 ];
ub = [5.15e-3  0.9     5     3e8   Tdata2(n2+1) Tdata3(n3+1)-5e6 Tdata4(n4+1) 1e7];
x  = [5.15e-3 0.006577733464 2.558775384581 1e8 ...
      Tdata2(n2)     Tdata3(n3)-6e6 Tdata4(n4)      1e6];

% RT timing (or dose-time parameter)
d1 = 0;
d2 = 4;
d3 = 8;
d4 = 15;

options = optimset('Disp','iter','TolX',1e-12,'TolFun',1e-12);

x = lsqcurvefit(@(x, td) cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2, d3, d4]), ...
                x, td, Tdata, lb, ub, options);

% Evaluate cost and metrics
cost_value = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2, d3, d4]);
Rd   = norm(cost_value - Tdata) / norm(Tdata);

% AIC / BIC
error      = cost_value - Tdata;
gaussian_norm = sqrt(sum(error.^2));    % same as norm(error)

Tdata_norm = norm(Tdata);
Rd2        = gaussian_norm / Tdata_norm;

n_total    = numel(Tdata);
num_params = numel(x);
AIC = n_total * log(sum(error.^2)/n_total) + 2 * num_params;
BIC = n_total * log(sum(error.^2)/n_total) + num_params * log(n_total);

names = {'\alpha','\beta','\lambda', ...
         'c_1(0)','c_2(0)','c_3(0)','c_4(0)','E_(0)'};
T = table(names(:), x(:), 'VariableNames', {'Parameter','Estimate'});
disp('================ FIT SUMMARY ================');
fprintf('AIC = %.3f | BIC = %.3f | Rd = %.4f | Rd2 = %.4f\n', AIC, BIC, Rd, Rd2);
disp(T);
t1=8;

save('../Results/fit_step3.mat','x');



% ---------- Local function ----------
function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, d)
    [~, Y1_common] = ode23s(@(t, y) raditation_model(t, y, x, d(1)), tdata(n1:m1), [x(4); x(8)]);
    [~, Y2_common] = ode23s(@(t, y) raditation_model(t, y, x, d(2)), tdata(n2:m2), [x(5); x(8)]);
    [~, Y3_common] = ode23s(@(t, y) raditation_model(t, y, x, d(3)), tdata(n3:m3), [x(6); x(8)]);
    [~, Y4_common] = ode23s(@(t, y) raditation_model(t, y, x, d(4)), tdata(n4:m4), [x(7); x(8)]);

    cost = [Y1_common(:,1) + Y1_common(:,2);
            Y2_common(:,1) + Y2_common(:,2);
            Y3_common(:,1) + Y3_common(:,2);
            Y4_common(:,1) + Y4_common(:,2)];

    assert(size(cost, 1) == size(tdata(n1:m1), 1) + size(tdata(n2:m2), 1) + ...
                           size(tdata(n3:m3), 1) + size(tdata(n4:m4), 1), ...
           'Inconsistent dimensions');
end
