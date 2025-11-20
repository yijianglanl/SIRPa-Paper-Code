
format long 

warning('off','all')
sensitivity_on = true;
% Upper bounds for datasets
m1 = 11;
m2 = 13;
m3 = 13;
m4=13;
% Lower bounds for datasets
n1 =5;
n2 = 5;
n3 = 5;
n4=5;
load('../Data/it2dose.txt');

tdata = it2dose(:, 1);
Tdata1 = it2dose(:, 2);
Tdata2 = it2dose(:, 3);
Tdata3 =it2dose(:, 4);
Tdata4 =it2dose(:, 5);
st1 = it2dose(:, 6);
st2 = it2dose(:, 7);
st3 =it2dose(:, 8);
st4 =it2dose(:, 9);

td = [tdata(n1:m1); tdata(n2:m2); tdata(n3:m3); tdata(n4:m4)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2); Tdata3(n3:m3); Tdata4(n4:m4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
d1 = 0;
d2 = 8;

lb = [1 2.3e-4 1 1 1 1 1e5 1e5 1 1 1 ];
ub = [30 2.3e-1 1e8 1e8 1e8 1e8 1e7 1e7 20 25 1e6 ];
x = [15 1e-2 1e6 1e6 1e6 1e6 1e6 1e6 5 19 1e4];

x = lsqcurvefit(@(x,td) cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2]), x, td, Tdata, lb, ub, ...
    optimset('Disp', 'iter', 'TolX', 10^(-15), 'TolFun', 10^(-15)));
cost_value = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2]);

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
% color options
color2 = [0.8, 0.2, 0.2]; 
color4 = [0.2, 0.6, 0.2];  
color3 = [0.2, 0.2, 0.6];  
color1 = [0.2, 0.2, 0.2];                    

 figure1 = figure;
set(figure1, 'Position', [100, 100, 800, 400]); 
% Plot
plot(tdata(n1:m1), cost_value(1:m1-n1+1)./1e6, '-', 'Color', color1, 'LineWidth', 3);
hold on;
errorbar(tdata(n1:m1), Tdata1(n1:m1)./1e6, st1(n1:m1)./1e6, 'o', 'Color', color1, 'MarkerSize', 8, 'MarkerFaceColor', color1, 'LineStyle', 'none', 'CapSize', 3, 'LineWidth', 2);
hold on;
plot(tdata(n2:m2), cost_value(m1-n1+2:m1-n1+m2-n2+2)./1e6, '-', 'Color', color2, 'LineWidth', 3);

errorbar(tdata(n2:m2), Tdata2(n2:m2)./1e6, st2(n2:m2)./1e6, 'o', 'Color', color2, 'MarkerSize', 8, 'MarkerFaceColor', color2, 'LineStyle', 'none', 'CapSize', 3, 'LineWidth', 2);
plot(tdata(n3:m3), cost_value(m1-n1+m2-n2+3:m1-n1+m2-n2+m3-n3+3)./1e6, '-', 'Color', color3, 'LineWidth', 3);

errorbar(tdata(n3:m3), Tdata3(n3:m3)./1e6, st3(n3:m3)./1e6, 'o', 'Color', color3, 'MarkerSize', 8, 'MarkerFaceColor', color3, 'LineStyle', 'none', 'CapSize', 3, 'LineWidth', 2);
plot(tdata(n4:m4), cost_value(m1-n1+m2-n2+m3-n3+4:m1-n1+m2-n2+m3-n3+m4-n4+4)./1e6, '-', 'Color', color4, 'LineWidth',3);

errorbar(tdata(n4:m4), Tdata4(n4:m4)./1e6, st4(n4:m4)./1e6, 'o', 'Color', color4, 'MarkerSize', 8, 'MarkerFaceColor', color4, 'LineStyle', 'none', 'CapSize', 3, 'LineWidth', 2);

% Adding annotations and a vertical line for treatment indication
xline(12, '--m', 'LineWidth', 2, 'Label', '8Gy Treatment','fontsize',19, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'middle');
xline(14, '--m', 'LineWidth', 2, 'Label', '8Gy Treatment','fontsize',19, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'middle');
grid on
% 
% Configuring axes labels and limits
xlim([7.5,24.3])
ylim([-0.11e8/1e6, 1.62e9/1e6]); % Adjust limits based on data
set(gca, 'FontSize', 19, 'LineWidth', 1.3);
%legend('Control data', 'Model Prediction ', 'WT 8Gy data  ','Model Prediction','8Gy + i.t BMDM data','Model Prediction' ,'8Gy+ i.t SIRPA BMDM data','Model Prediction' ,'Location', 'westoutside', 'FontSize', 20);

xlabel('\fontsize{20}Time (days)',FontWeight='bold');
ylabel('\fontsize{20}Tumor volume (mm^3)',FontWeight='bold');

function dydt = g1(t, y, x, d)
    dydt = zeros(2,1);
    t0 = d/1.2;
    if (d == 0)
        u = 1;
    else
          u=2*(2.558775384581*t0+ exp(-2.558775384581*t0)- 1)/((2.558775384581^2)*t0*t0);
    end
      
       S1=exp(-0.005147638197*d - 0.006577733464*d*d*u);
    
    dydt(1) = 0.453709396815* y(1) * (1 - y(1) /1759812730.672900915146) -0.000000044133* y(1) * y(2);
    dydt(2) =16338.485039659030 -   0.069125930809* y(2) - 0.000000001991* y(1) * y(2);
   
    if (t >= 12 && t < 14)

        dydt(1) = dydt(1)*S1 ;
        dydt(2) = dydt(2);
        
    elseif (t>=14)
        dydt(1) = dydt(1)*S1;
        dydt(2) = dydt(2);
        
    else
        dydt(1) = dydt(1);
        dydt(2) = dydt(2);
        
    end
end
function dydt = g2(t, y, x, d)
    dydt = zeros(3,1);
    t0 = d/1.2;
    if (d == 0)
        u = 1;
    else
         u=2*(2.558775384581*t0+ exp(-2.558775384581*t0)- 1)/((2.558775384581^2)*t0*t0);
    end
      
       S1=exp(-0.005147638197*d - 0.006577733464*d*d*u);
   
    dydt(1) = 0.453709396815* y(1) * (1 - y(1) /1759812730.672900915146) -0.000000044133* y(1) * y(2)-4.43e-8*y(1)*y(3);
    dydt(2) =16338.485039659030 -   0.069125930809* y(2) - 0.000000001991* y(1) * y(2);
    dydt(3) = 0.000000439623   -0.394943790649 *y(3)-0.000000815338*y(1)*y(3);

   if (t >= 12 && t < 14)

        dydt(1) = dydt(1)*S1 - dydt(1)*(1-S1)*8.2474254899;
        dydt(2) = dydt(2);
        dydt(3) = dydt(3);
    elseif (t>=14)
        dydt(1) = dydt(1)- dydt(1)*(1-S1)*19.9;
        dydt(2) = dydt(2);
        dydt(3) = dydt(3);
    else
        dydt(1) = dydt(1);
        dydt(2) = dydt(2);
        dydt(3) = dydt(3);
    end
end

function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, d)
    [~, Y1_common] = ode23s(@(t, y) g1(t, y, x, d(1)), tdata(n1:m1), [x(3); x(7)]);
    [~, Y2_common] = ode23s(@(t, y) g1(t, y, x, d(2)), tdata(n2:m2), [x(4); x(7)]);
    [~, Y3_common] = ode23s(@(t, y) g2(t, y, x, d(1)), tdata(n3:m3), [x(5); x(11); x(8)]);
    [~, Y4_common] = ode23s(@(t, y) g2(t, y, x, d(2)), tdata(n4:m4), [x(6); x(11); x(8)]);

    % Concatenating the sums into a single vector
    cost = [Y1_common(:, 1) + Y1_common(:, 2) ;
            Y2_common(:, 1) + Y2_common(:, 2) ;
            Y3_common(:, 1) + Y3_common(:, 2) + Y3_common(:, 3);
            Y4_common(:, 1) + Y4_common(:, 2) + Y4_common(:, 3)];

    % size of the cost matches the size of Tdata
    assert(size(cost, 1) == size(tdata(n1:m1), 1) + size(tdata(n2:m2), 1) + size(tdata(n3:m3), 1) + size(tdata(n4:m4), 1), 'Inconsistent dimensions');
end
