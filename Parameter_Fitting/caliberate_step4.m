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

% ===== Calibrate step 4 =====
m1 = 10; m2 = 12; m3 = 13; m4 = 13;
n1 = 5;  n2 = 4;  n3 = 4;  n4 = 4;

load('../Data/MCSrd100_4.txt');

tdata  = MCSrd100_4(:, 1);
Tdata1 = MCSrd100_4(:, 2);
Tdata2 = MCSrd100_4(:, 3);
Tdata3 = MCSrd100_4(:, 4);
Tdata4 = MCSrd100_4(:, 5);

st1 = MCSrd100_4(:, 6);
st2 = MCSrd100_4(:, 7);
st3 = MCSrd100_4(:, 8);
st4 = MCSrd100_4(:, 9);

lb = [0 0 0    1000 1000 1000 1000 1e3 1e3];
ub = [1000 1000 1000  Tdata1(n1+1) Tdata2(n2+1) Tdata3(n3+1) Tdata4(n4+1) 1e7 1e7];
x  = [5 5 5 Tdata1(n1) Tdata2(n2) Tdata3(n3) Tdata4(n4) 1e6 1e6];
d1 = 0; d2 = 4; d3 = 8; d4 = 15;
td    = [tdata(n1:m1); tdata(n2:m2); tdata(n3:m3); tdata(n4:m4)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2); Tdata3(n3:m3); Tdata4(n4:m4)];
options = optimset('Disp', 'iter', 'TolX', 1e-12, 'TolFun', 1e-12);

x = lsqcurvefit(@(x, td) cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2, d3, d4]), ...
                x, td, Tdata, lb, ub, options);

cost_value = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2, d3, d4]);
Rd = norm(cost_value - Tdata) / norm(Tdata);

residuals  = cost_value - Tdata;
n_total    = numel(Tdata);
num_params = numel(x);
AIC = n_total*log(sum(residuals.^2)/n_total) + 2*num_params;
BIC = n_total*log(sum(residuals.^2)/n_total) + num_params*log(n_total);

names = {'A_4','A_8','A_{15}', ...
         'c_1(0)','c_2(0)','c_3(0)','c_4(0)','E_(0)','M(0)'};
T = table(names(:), x(:), 'VariableNames', {'Parameter','Estimate'});
disp('================ FIT SUMMARY ================');
fprintf('AIC = %.3f | BIC = %.3f | Rd = %.4f\n', AIC, BIC, Rd);
disp(T);
save('../Results/fit_step4.mat','x');


function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, d)
    [~, Y1_common] = ode23s(@(t, y) ICD_model(t, y, x,d(1)), tdata(n1:m1), [x(4); x(8); x(9)]);
    [~, Y2_common] = ode23s(@(t, y) ICD_model(t, y, x, d(2)), tdata(n2:m2), [x(5); x(8); x(9)]);
    [~, Y3_common] = ode23s(@(t, y) ICD_model(t, y, x, d(3)), tdata(n3:m3), [x(6); x(8); x(9)]);
    [~, Y4_common] = ode23s(@(t, y) ICD_model(t, y, x, d(4)), tdata(n4:m4), [x(7); x(8); x(9)]);

   
    cost = [Y1_common(:, 1) + Y1_common(:, 2) + Y1_common(:, 3);
            Y2_common(:, 1) + Y2_common(:, 2) + Y2_common(:, 3);
            Y3_common(:, 1) + Y3_common(:, 2) + Y3_common(:, 3);
            Y4_common(:, 1) + Y4_common(:, 2) + Y4_common(:, 3)];

    assert(size(cost, 1) == size(tdata(n1:m1), 1) + size(tdata(n2:m2), 1) + size(tdata(n3:m3), 1) + size(tdata(n4:m4), 1), 'Inconsistent dimensions');
end
