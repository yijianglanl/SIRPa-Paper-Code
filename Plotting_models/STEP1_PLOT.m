% ====== Plot Step 1 (no fitting, using precomputed x) ======
clear; clc;
format long
warning('off','all')
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
% ---------------- Index windows ----------------
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


S = load('../Results/WT_fit_step1.mat');  % <-- adjust path/name
x = S.x;   
% ---------------- Compute model output with fixed x ----------------
yfit = cost_function(x,tdata,n1,n2,n3,n4,m1,m2,m3,m4);

% ---------------- (Optional) Fit metrics for reporting ----------------
Rd   = norm(yfit - Tdata)/norm(Tdata);
residuals  = yfit - Tdata;
n_total    = numel(Tdata);
num_params = numel(x);
AIC = n_total*log(sum(residuals.^2)/n_total) + 2*num_params;
BIC = n_total*log(sum(residuals.^2)/n_total) + num_params*log(n_total);

fprintf('AIC: %.3f\nBIC: %.3f\nRd: %.4f\n', AIC, BIC, Rd);
fprintf('Params (x1..x11):\n'); disp(x(:).')

names = {'c_1','c_{max}','\phi_e','\gamma_e','\delta_e','\eta_e', ...
         'c_1(0)','c_2(0)','c_3(0)','c_4(0)','E_(0)'};
T = table(names(:), x(:), 'VariableNames', {'Parameter','Estimate'});
disp('================ FIT SUMMARY ================');
fprintf('AIC = %.3f | BIC = %.3f | Rd = %.4f\n', AIC, BIC, Rd);
disp(T);

% ---------------- Plot: data vs model (scaled by 1e6) ----------------
L1 = m1-n1+1; L2 = m2-n2+1; L3 = m3-n3+1; L4 = m4-n4+1;
i1 = 1:L1;
i2 = L1 + (1:L2);
i3 = L1+L2 + (1:L3);
i4 = L1+L2+L3 + (1:L4);

color2 = [0.8, 0.2, 0.2];  
color4 = [0.2, 0.6, 0.2];  
color3 = [0.2, 0.2, 0.6]; 
color1 = [0.2, 0.2, 0.2];  

figure; hold on; box on;

h1 = plot(tdata(n1:m1), Tdata1(n1:m1)/1e6,  '<', 'Color', color1, 'LineStyle','none', 'MarkerSize',6, 'LineWidth',1.5);
h2 = plot(tdata(n1:m1), yfit(i1)/1e6, '-', 'Color', color1, 'LineWidth',2.5);

h3 = plot(tdata(n2:m2), Tdata2(n2:m2)/1e6,  'o', 'Color', color2, 'LineStyle','none', 'MarkerSize',6, 'LineWidth',1.5);
h4 = plot(tdata(n2:m2), yfit(i2)/1e6, '-', 'Color', color2, 'LineWidth',2.5);

h5 = plot(tdata(n3:m3), Tdata3(n3:m3)/1e6,  '*', 'Color', color3, 'LineStyle','none', 'MarkerSize',6, 'LineWidth',1.5);
h6 = plot(tdata(n3:m3), yfit(i3)/1e6, '-', 'Color', color3, 'LineWidth',2.5);

h7 = plot(tdata(n4:m4), Tdata4(n4:m4)/1e6,  's', 'Color', color4, 'LineStyle','none', 'MarkerSize',6, 'LineWidth',1.5);
h8 = plot(tdata(n4:m4), yfit(i4)/1e6, '-', 'Color', color4, 'LineWidth',2.5);

xlabel('Time (days)','FontSize',13);
ylabel('Tumor Volume (10^6 mm^3)','FontSize',13);
xlim([min([tdata(n1),tdata(n2),tdata(n3),tdata(n4)])-0.5, ...
      max([tdata(m1),tdata(m2),tdata(m3),tdata(m4)])+0.5]);
set(gca,'FontSize',12);
legend([h1 h2 h3 h4 h5 h6 h7 h8], ...
       {'Data set 1','Model fit 1','Data set 2','Model fit 2', ...
        'Data set 3','Model fit 3','Data set 4','Model fit 4'}, ...
       'Location','northwest','FontSize',10);
grid on; hold off;

% ---------------- Local ODE and cost function ----------------
function dydt = tumor_growth(~,y,x)
    dydt = zeros(2,1);
    dydt(1) = x(1)*y(1)*(1 - y(1)/x(2)) - x(3)*y(1)*y(2);
    dydt(2) = x(4) - x(5)*y(2) - x(6)*y(1)*y(2);
end

function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4)
    [~,Y1] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n1:m1), [x(7);  x(11)]);
    [~,Y2] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n2:m2), [x(8);  x(11)]);
    [~,Y3] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n3:m3), [x(9);  x(11)]);
    [~,Y4] = ode23s(@(t,y) tumor_growth(t,y,x), tdata(n4:m4), [x(10); x(11)]);

    Y1s = Y1(:,1);
    Y2s = Y2(:,1);
    Y3s = Y3(:,1);
    Y4s = Y4(:,1);

    cost = [Y1s; Y2s; Y3s; Y4s];
end
