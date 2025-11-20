% ====== Calibrate step 2  ======
format long
clear; clc;
warning('off','all')

clearvars -global


mf = mfilename('fullpath');
if isempty(mf)
    
    project_root = pwd;
else
    project_root = fileparts(mf);
end


paths.models    = fullfile(project_root, 'Models');
addpath(genpath(paths.models));

%% update 
addpath('../Models'); %% adjust for your computer

m1 = 17;
m2 = 14;
m3 = 11;
m4 = 8;

% Lower bounds for datasets

n1 = 10;
n2 = 7;
n3 = 4;
n4 = 2;

load('../Data/MCS_cells.txt'); %% adjust for your computer

tdata = MCS_cells(:,1);
Tdata1 = MCS_cells(:,2);
Tdata2 = MCS_cells(:,3);
Tdata3 = MCS_cells(:,4);
Tdata4 = MCS_cells(:,5);

td    = [tdata(n1:m1); tdata(n2:m2); tdata(n3:m3); tdata(n4:m4)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2); Tdata3(n3:m3); Tdata4(n4:m4)];

lb = [1e-7  1e-8 2e-3 1e-8 1000 1000 1000 1000 1 1];

ub = [1e-5   1e-6 8e-1 1e-6 5e6 1e8 1e8 1e8 1e6 1e6];

x  = [9e-7   4e-7 3e-1 8e-7 1e5 1e6 1e6 1e6 1e4 1e4];

x = lsqcurvefit(@(x,td) cost_function(x,tdata,n1,n2,n3,n4,m1,m2,m3,m4),x,td,Tdata,lb,ub,...
                             optimset('Disp','iter','TolX',10^(-15),'TolFun',10^(-15)));
                         

% ---------------- Fit metrics ----------------
yfit = cost_function(x,tdata,n1,n2,n3,n4,m1,m2,m3,m4);
Rd   = norm(yfit - Tdata)/norm(Tdata);

residuals  = yfit - Tdata;
n_total    = numel(Tdata);
num_params = numel(x);
AIC = n_total*log(sum(residuals.^2)/n_total) + 2*num_params;
BIC = n_total*log(sum(residuals.^2)/n_total) + num_params*log(n_total);

fprintf('AIC: %.3f\nBIC: %.3f\nRd: %.4f\n', AIC, BIC, Rd);
disp('Final parameters x:');
disp(x(:).')

names = {'\phi_m','\gamma_m','\delta_m','\eta_m', ...
         'c_1(0)','c_2(0)','c_3(0)','c_4(0)','E_(0)','M_(0)'};
T = table(names(:), x(:), 'VariableNames', {'Parameter','Estimate'});
disp('================ FIT SUMMARY ================');
fprintf('AIC = %.3f | BIC = %.3f | Rd = %.4f\n', AIC, BIC, Rd);
disp(T);
save('../Results/fit_step2.mat','x');

% ====== Local function ======
function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4)
    [~, Y1_common] = ode23s(@(t, y) MC38_SIRPA_KO(t, y, x), tdata(n1:m1), [x(5); x(9); x(10)]);
    [~, Y2_common] = ode23s(@(t, y) MC38_SIRPA_KO(t, y, x), tdata(n2:m2), [x(6);x(9); x(10)]);
    [~, Y3_common] = ode23s(@(t, y) MC38_SIRPA_KO(t, y, x), tdata(n3:m3), [x(7); x(9); x(10)]);
    [~, Y4_common] = ode23s(@(t, y) MC38_SIRPA_KO(t, y, x), tdata(n4:m4), [x(8); x(9); x(10)]);

    cost = [Y1_common(:, 1)+Y1_common(:, 2)+Y1_common(:, 3);
            Y2_common(:, 1)+Y2_common(:, 2)+Y2_common(:, 3);
            Y3_common(:, 1)+Y3_common(:, 2)+Y3_common(:, 3);
            Y4_common(:, 1)+Y4_common(:, 2)+Y4_common(:, 3)];
end
