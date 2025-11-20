format long
warning('off','all')
clear; clc; close all
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

td    = [tdata(n1:m1); tdata(n2:m2); tdata(n3:m3); tdata(n4:m4)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2); Tdata3(n3:m3); Tdata4(n4:m4)];

lb = [0 0 0    1000 1000 1000 1000 1e3 1e3];
ub = [1000 1000 1000  Tdata1(n1+1) Tdata2(n2+1) Tdata3(n3+1) Tdata4(n4+1) 1e7 1e7];
x0 = [5 5 5 Tdata1(n1) Tdata2(n2) Tdata3(n3) Tdata4(n4) 1e6 1e6];

d1 = 0; d2 = 4; d3 = 8; d4 = 15;

S = load('../Results/fit_step4.mat');
x = S.x;

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

t1 = 8;

L1 = m1 - n1 + 1;
L2 = m2 - n2 + 1;
L3 = m3 - n3 + 1;
L4 = m4 - n4 + 1;

i1 = 1:L1;
i2 = L1 + (1:L2);
i3 = L1+L2 + (1:L3);
i4 = L1+L2+L3 + (1:L4);

c2 = [0.8, 0.2, 0.2];
c4 = [0.2, 0.6, 0.2];
c3 = [0.2, 0.2, 0.6];
c1 = [0.2, 0.2, 0.2];

figure; hold on; box on;

h1 = plot(tdata(n1:m1), Tdata1(n1:m1)/1e6, 'o', ...
          'Color', c1, 'MarkerFaceColor', c1, ...
          'LineStyle','none', 'MarkerSize', 5);
h2 = plot(tdata(n1:m1), cost_value(i1)/1e6, '-', ...
          'Color', c1, 'LineWidth', 2.2);
if exist('st1','var')
    errorbar(tdata(n1:m1), Tdata1(n1:m1)/1e6, st1(n1:m1)/1e6, ...
             'Color', c1, 'LineStyle','none', 'CapSize', 3, 'LineWidth', 1.2);
end

h3 = plot(tdata(n2:m2), Tdata2(n2:m2)/1e6, 's', ...
          'Color', c2, 'MarkerFaceColor', c2, ...
          'LineStyle','none', 'MarkerSize', 5);
h4 = plot(tdata(n2:m2), cost_value(i2)/1e6, '-', ...
          'Color', c2, 'LineWidth', 2.2);
if exist('st2','var')
    errorbar(tdata(n2:m2), Tdata2(n2:m2)/1e6, st2(n2:m2)/1e6, ...
             'Color', c2, 'LineStyle','none', 'CapSize', 3, 'LineWidth', 1.2);
end

h5 = plot(tdata(n3:m3), Tdata3(n3:m3)/1e6, 'd', ...
          'Color', c3, 'MarkerFaceColor', c3, ...
          'LineStyle','none', 'MarkerSize', 5);
h6 = plot(tdata(n3:m3), cost_value(i3)/1e6, '-', ...
          'Color', c3, 'LineWidth', 2.2);
if exist('st3','var')
    errorbar(tdata(n3:m3), Tdata3(n3:m3)/1e6, st3(n3:m3)/1e6, ...
             'Color', c3, 'LineStyle','none', 'CapSize', 3, 'LineWidth', 1.2);
end

h7 = plot(tdata(n4:m4), Tdata4(n4:m4)/1e6, '^', ...
          'Color', c4, 'MarkerFaceColor', c4, ...
          'LineStyle','none', 'MarkerSize', 5);
h8 = plot(tdata(n4:m4), cost_value(i4)/1e6, '-', ...
          'Color', c4, 'LineWidth', 2.2);
if exist('st4','var')
    errorbar(tdata(n4:m4), Tdata4(n4:m4)/1e6, st4(n4:m4)/1e6, ...
             'Color', c4, 'LineStyle','none', 'CapSize', 3, 'LineWidth', 1.2);
end

xlabel('Time (days)','FontSize',13);
ylabel('Tumor Volume (10^6 mm^3)','FontSize',13);

xlim([min([tdata(n1),tdata(n2),tdata(n3),tdata(n4)])-0.5, ...
      max([tdata(m1),tdata(m2),tdata(m3),tdata(m4)])+0.5]);

set(gca,'FontSize',12);
legend([h1 h2 h3 h4 h5 h6 h7 h8], ...
       {'0 Gy data','0 Gy model','4 Gy data','4 Gy model', ...
        '8 Gy data','8 Gy model','15 Gy data','15 Gy model'}, ...
       'Location','northwest','FontSize',10);
grid on; hold off;

function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, d)
    [~, Y1_common] = ode23s(@(t, y) ICD_model(t, y, x, d(1)), tdata(n1:m1), [x(4); x(8); x(9)]);
    [~, Y2_common] = ode23s(@(t, y) ICD_model(t, y, x, d(2)), tdata(n2:m2), [x(5); x(8); x(9)]);
    [~, Y3_common] = ode23s(@(t, y) ICD_model(t, y, x, d(3)), tdata(n3:m3), [x(6); x(8); x(9)]);
    [~, Y4_common] = ode23s(@(t, y) ICD_model(t, y, x, d(4)), tdata(n4:m4), [x(7); x(8); x(9)]);

    cost = [Y1_common(:, 1) + Y1_common(:, 2) + Y1_common(:, 3);
            Y2_common(:, 1) + Y2_common(:, 2) + Y2_common(:, 3);
            Y3_common(:, 1) + Y3_common(:, 2) + Y3_common(:, 3);
            Y4_common(:, 1) + Y4_common(:, 2) + Y4_common(:, 3)];
end
