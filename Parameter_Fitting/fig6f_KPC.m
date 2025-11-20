
format long 

warning('off','all')
load('../Data/KPCRt.txt');

m1 = 14;
m2 = 15;
m3 = 14;
m4=15;
% Lower bounds for datasets
n1 =8;
n2 = 8;
n3 = 8;
n4=7;

tdata =  KPCRt(:, 1);
Tdata1 = KPCRt(:, 2);
Tdata2 = KPCRt(:, 3);
Tdata3 = KPCRt(:, 4);
Tdata4 = KPCRt(:, 5);
st1 = KPCRt(:, 6);
st2 =KPCRt(:, 7);
st3 =KPCRt(:, 8);
st4 =KPCRt(:, 9);
td = [tdata(n1:m1); tdata(n2:m2); tdata(n3:m3); tdata(n4:m4)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2); Tdata3(n3:m3); Tdata4(n4:m4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
d1 = 0;
d2 = 8;

lb = [0.3 8e8 1 1 1 1 1e5 1e5 0 0];
ub = [0.65 5e9 1e8 1e8 1e8 1e8 1e7 1e7 100 1e6];
x = [0.45 2e9 1e6 1e6 1e6 1e6 1e6 1e6 5 1e5];



x = lsqcurvefit(@(x,td) cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2]), x, td, Tdata, lb, ub, ...
    optimset('Disp', 'iter', 'TolX', 10^(-15), 'TolFun', 10^(-15)));
cost_value = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2]);
 %cost_value = cost_function(x, tdata, n1, n2, n3,n4, m1, m2, m3,m4,[d1, d2, d3,d4]);
    Rd = norm(cost_value - Tdata) / norm(Tdata);
residuals = cost_value - Tdata;
n_total = numel(Tdata);
num_params = numel(x);
 disp('____________________________________________________________________________________________________')
    disp('..alpha..Beta...lambda...S2...c1...c2...c3...c4..d1...d2...d3...d4...m0. Rd..')
    disp('____________________________________________________________________________________________________')
    
    % CREATION OF THE VECTOR FOR OUTPUT
    fprintf('%10.10f %10.10f %10.5f %10.10f %10.10f %10.5f %10.10f %10.10f %10.10f %10.5f %10.10f %10.5f %10.5f %10.5f\n', ...
        x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8),x(9),x(10),Rd);
    

%--------------------------------------------------------------------
%  color options
color2 = [0.8, 0.2, 0.2];  
color4 = [0.2, 0.6, 0.2];  
color3 = [0.2, 0.2, 0.6];  
color1 = [0.2, 0.2, 0.2];                   
figure1 = figure;
set(figure1, 'Position', [100, 100, 800, 400]); 

% Ploting all datasets and model results 
errorbar(tdata(n1:m1), Tdata1(n1:m1)./1e6, st1(m1:n1)./1e6, 'o', 'Color', color1, 'MarkerSize', 8, 'MarkerFaceColor', color1, 'LineStyle', 'none', 'CapSize', 3, 'LineWidth', 2);
hold on;
plot(tdata(n1:m1), cost_value(1:m1-n1+1)./1e6, '-', 'Color', color1, 'LineWidth', 3);

errorbar(tdata(n2:m2), Tdata2(n2:m2)./1e6, st2(n2:m2)./1e6, 'o', 'Color', color2, 'MarkerSize', 8, 'MarkerFaceColor', color2, 'LineStyle', 'none', 'CapSize', 3, 'LineWidth', 2);
plot(tdata(n2:m2), cost_value(m1-n1+2:m1-n1+m2-n2+2)./1e6, '-', 'Color', color2, 'LineWidth', 3);

errorbar(tdata(n3:m3), Tdata3(n3:m3)./1e6, st3(n3:m3)./1e6, 'o', 'Color', color3, 'MarkerSize', 8, 'MarkerFaceColor', color3, 'LineStyle', 'none', 'CapSize', 3, 'LineWidth', 2);
plot(tdata(n3:m3), cost_value(m1-n1+m2-n2+3:m1-n1+m2-n2+m3-n3+3)./1e6, '-', 'Color', color3, 'LineWidth', 3);

errorbar(tdata(n4:m4), Tdata4(n4:m4)./1e6, st4(n4:m4)./1e6, 'o', 'Color', color4, 'MarkerSize', 8, 'MarkerFaceColor', color4, 'LineStyle', 'none', 'CapSize', 3, 'LineWidth', 2);
plot(tdata(n4:m4), cost_value(m1-n1+m2-n2+m3-n3+4:m1-n1+m2-n2+m3-n3+m4-n4+4)./1e6, '-', 'Color', color4, 'LineWidth', 3);


ylim([0.1e8./1e6, 1.75e9./1e6]);
grid on;
% Adding annotations and a vertical line for treatment indication
xline(18, ':m', 'LineWidth', 2.5, 'Label', '8Gy Treatment','fontsize',17, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'middle');
%xline(14, ':m', 'LineWidth', 2.5, 'Label', 'Treatment Day','fontsize',15, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'middle');

% Configuring axes labels and limits
xlabel('Time (days)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Tumor volume (mm^3)', 'FontSize', 18, 'FontWeight', 'bold');

xlim([11.5,28.5]);
set(gca, 'FontSize', 19, 'LineWidth', 1.3);
grid on
%legend({'Control data', 'Model Fit', 'WT+8Gy data', 'Model Fit', 'SIRP\alpha-/- 0Gy data', 'Model Fit', 'SIRP\alpha-/- 8Gy data', 'Model Fit'},...
       % 'FontSize', 23, 'Location', 'northwest');

%annotation('textbox', [0.05, 0.85, 0.1, 0.1], 'String', 'B', 'EdgeColor', 'none', 'FontSize', 17, 'FontWeight', 'bold');



function dydt = g(t, y, x)
    dydt = zeros(2,1);
    
    dydt(1) = x(1)* y(1) * (1 - y(1) /x(2)) -0.000000044133* y(1) * y(2);
    dydt(2) =16338.485039659030 -   0.069125930809* y(2) - 0.000000001991* y(1) * y(2);
 
end
function dydt = g1(t, y, x, d)
    dydt = zeros(2,1);
    t0 = d/1.2;
    if (d == 0)
        u = 1;
    else
        u=2*(2.558775384581*t0+ exp(-2.558775384581*t0)- 1)/((2.558775384581^2)*t0*t0);
    end
      
       S1=exp(-0.005147638197*d - 0.006577733464*d*d*u);
    dydt(1) = x(1)* y(1) * (1 - y(1) /x(2)) -0.000000044133* y(1) * y(2);
    dydt(2) =16338.485039659030 -   0.069125930809* y(2) - 0.000000001991* y(1) * y(2);

    if (t >= 18)
        dydt(1) = dydt(1)*S1 ;
        dydt(2) = dydt(2);
        
    else
        dydt(1) = dydt(1);
        dydt(2) = dydt(2);
       
    end
end
function dydt = g2(t, y, x)
    dydt = zeros(3,1);
   
   
    dydt(1) = x(1)* y(1) * (1 - y(1) /x(2)) -0.000000044133* y(1) * y(2)-0.000000907400*y(1)*y(3);
    dydt(2) =16338.485039659030 -   0.069125930809* y(2) - 0.000000001991* y(1) * y(2);
    dydt(3) = 0.000000439623   -0.394943790649 *y(3)-0.000000815338*y(1)*y(3);

 
end
function dydt = g3(t, y, x, d)
    dydt = zeros(3,1);
    t0 = d/1.2;
    if (d == 0)
        u = 1;
    else
        u=2*(2.558775384581*t0+ exp(-2.558775384581*t0)- 1)/((2.558775384581^2)*t0*t0);
    end
      
       S1=exp(-0.005147638197*d - 0.006577733464*d*d*u);
    A = x(9);
    dydt(1) = x(1)* y(1) * (1 - y(1) /x(2)) -0.000000044133* y(1) * y(2)-0.000000907400*y(1)*y(3);
    dydt(2) =16338.485039659030 -   0.069125930809* y(2) - 0.000000001991* y(1) * y(2);
    dydt(3) = 0.000000439623   -0.394943790649 *y(3)-0.000000815338*y(1)*y(3);

   if (t >= 18 )

        dydt(1) = dydt(1)*S1 - dydt(1)*(1-S1)*A;
        dydt(2) = dydt(2);
        dydt(3) = dydt(3);
    
    else
        dydt(1) = dydt(1);
        dydt(2) = dydt(2);
        dydt(3) = dydt(3);
    end
end
function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, d)
    [~, Y1_common] = ode23s(@(t, y) g(t, y, x), tdata(n1:m1), [x(3); x(7)]);
    [~, Y2_common] = ode23s(@(t, y) g1(t, y, x, d(2)), tdata(n2:m2), [x(4); x(7)]);
    [~, Y3_common] = ode23s(@(t, y) g2(t, y, x), tdata(n3:m3), [x(5); x(10); x(8)]);
    [~, Y4_common] = ode23s(@(t, y) g3(t, y, x, d(2)), tdata(n4:m4), [x(6); x(10); x(8)]);

    cost = [Y1_common(:, 1) + Y1_common(:, 2) ;
            Y2_common(:, 1) + Y2_common(:, 2) ;
            Y3_common(:, 1) + Y3_common(:, 2) + Y3_common(:, 3);
            Y4_common(:, 1) + Y4_common(:, 2) + Y4_common(:, 3)];

    assert(size(cost, 1) == size(tdata(n1:m1), 1) + size(tdata(n2:m2), 1) + size(tdata(n3:m3), 1) + size(tdata(n4:m4), 1), 'Inconsistent dimensions');
end
