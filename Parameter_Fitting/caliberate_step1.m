% ====== Calibrate step 1 ======
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
load('../Data/WT.txt');

m1 = 16;
m2 = 14;
m3 = 11;
m4 = 8;

% Lower bounds for datasets

n1 = 8;
n2 = 7;
n3 = 4;
n4 = 2;

tdata  = WT(:,1);
Tdata1 = WT(:,2);
Tdata2 = WT(:,3);
Tdata3 = WT(:,4);
Tdata4 = WT(:,5);

td    = [tdata(n1:m1); tdata(n2:m2); tdata(n3:m3); tdata(n4:m4)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2); Tdata3(n3:m3); Tdata4(n4:m4)];

% ---------------- Bounds, initial guess ----------------
lb = [0.41 1.2e9     1e-8   1.3e3  1e-3    3.34e-11 10000 10000 10000 10000 10000 ];
ub = [0.45  2.6e9   1e-6  1.3e5  5e-1    3.34e-9  Tdata1(n1+1) Tdata2(n2+1) Tdata3(n3+1) Tdata4(n4+1) 1e7 ];
x0 = [0.42   2e9    1e-7   1.3e4  2e-2    3.34e-10   Tdata1(n1) Tdata2(n2) Tdata3(n3) Tdata4(n4) 1e6 ];


opts = optimset('Display','iter','TolX',1e-15,'TolFun',1e-15);

x = lsqcurvefit(@(x,td) cost_function(x,tdata,n1,n2,n3,n4,m1,m2,m3,m4), ...
                x0, td, Tdata, lb, ub, opts);

% ---------------- Fit metrics ----------------
yfit = cost_function(x,tdata,n1,n2,n3,n4,m1,m2,m3,m4);
Rd   = norm(yfit - Tdata)/norm(Tdata);

residuals  = yfit - Tdata;
n_total    = numel(Tdata);
num_params = numel(x);
AIC = n_total*log(sum(residuals.^2)/n_total) + 2*num_params;
BIC = n_total*log(sum(residuals.^2)/n_total) + num_params*log(n_total);

fprintf('AIC: %.3f\nBIC: %.3f\nRd: %.4f\n', AIC, BIC, Rd);
fprintf('Params (x1..x10, y1(0)..y4(0), c15, c16):\n');
disp(x(:).')

names = {'c_1','c_{max}','\phi_e','\gamma_e','\delta_e','\eta_e', ...
         'c_1(0)','c_2(0)','c_3(0)','c_4(0)','E_(0)'};
T = table(names(:), x(:), 'VariableNames', {'Parameter','Estimate'});
disp('================ FIT SUMMARY ================');
fprintf('AIC = %.3f | BIC = %.3f | Rd = %.4f\n', AIC, BIC, Rd);
disp(T);
%save('WT_fit_step1.mat', 'x');   % Save fitted params
save('../Results/WT_fit_step1.mat','x');


% ====== Local function using external tumor_growth.m ======
function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4)
    

    [~,Y1] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n1:m1), [x(7); x(11)]);
    [~,Y2] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n2:m2), [x(8); x(11)]);
    [~,Y3] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n3:m3), [x(9); x(11)]);
    [~,Y4] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n4:m4), [x(10); x(11)]);

    Y1s = Y1(:,1);  
    Y2s = Y2(:,1);
    Y3s = Y3(:,1);
    Y4s = Y4(:,1);

    cost = [Y1s; Y2s; Y3s; Y4s];
end
