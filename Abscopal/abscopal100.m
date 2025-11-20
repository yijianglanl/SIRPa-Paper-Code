
format long 

warning('off','all')

% Upper bounds for datasets
m1 = 10;
m2 = 16;
m3 = 16;
m4=10;

% Lower bounds for datasets
n1 = 5;
n2 = 5;
n3 = 6;
n4 = 6;

load('../Data/abs100.txt');
tdata =  abs100(:, 1);
Tdata1 = abs100(:, 2);
Tdata2 = abs100(:, 3);
Tdata3 = abs100(:, 4);
Tdata4 = abs100(:, 5);
st1 = abs100(:, 6);
st2 = abs100(:, 7);
st3 = abs100(:, 8);
st4 = abs100(:, 9);

td = [tdata(n1:m1); tdata(n2:m2); tdata(n3:m3); tdata(n4:m4)];
Tdata = [Tdata1(n1:m1); Tdata2(n2:m2); Tdata3(n3:m3); Tdata4(n4:m4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
d1 = 0;
d2 = 8;

lb = [1, 1, 1, 1, 1, 1,0.2,1, 0.1];
ub = [Tdata1(n1+1), Tdata2(n2+1), Tdata3(n3+1), Tdata4(n4+1), 1e7, 1e7,0.8,1e7,0.6];
x = [Tdata1(n1), Tdata2(n2), Tdata3(n3), Tdata4(n4), 1e6, 1e6,0.4,1e6,0.3];

% Fitting the model using lsqcurvefit
x = lsqcurvefit(@(x,td) cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2]), x, td, Tdata, lb, ub, ...
    optimset('Disp', 'iter', 'TolX', 10^(-15), 'TolFun', 10^(-15)));


cost_value = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, [d1, d2]);
Rd = norm(cost_value - Tdata) / norm(Tdata);


fprintf('%10.10f %10.10f %10.5f %10.10f %10.10f %10.5f %10.10f %10.10f %10.10f %10.5f %10.10f\n', ...
    x(1), x(2), x(3), x(4), x(5), x(6),  Rd);


color2 = [0.8, 0.2, 0.2];  % Dark Red
color4 = [0.2, 0.6, 0.2];  % Dark Green
color3 = [0.2, 0.2, 0.6];  % Dark Blue
color1 = [0.2, 0.2, 0.2];  % Dark Gray                   


figure1 = figure;

errorbar(tdata(n1:m1), Tdata1(n1:m1)./1e6, st1(n1:m1)./1e6, '-o', 'Color', color1, 'MarkerSize', 8, 'MarkerFaceColor', color1, 'LineStyle', '--', 'CapSize', 3, 'LineWidth', 1.5);
hold on;
plot(tdata(n2:m2), cost_value(m1-n1+2:m1-n1+m2-n2+2)./1e6, '-', 'Color', color2, 'LineWidth', 3);


errorbar(tdata(n2:m2), Tdata2(n2:m2)./1e6, st2(n2:m2)./1e6, 'o', 'Color', color2, 'MarkerSize', 8, 'MarkerFaceColor', color2, 'LineStyle', '--', 'CapSize', 3, 'LineWidth', 2);
%plot(tdata(n4:m4), , '-', 'Color', color4, 'LineWidth', 3);

errorbar(tdata(n4:m4), Tdata4(n4:m4)./1e6, st4(n4:m4)./1e6, '-o', 'Color', color4, 'MarkerSize', 8, 'MarkerFaceColor', color4, 'LineStyle', '--', 'CapSize', 3, 'LineWidth', 2);
plot(tdata(n3:m3), cost_value(m1-n1+m2-n2+3:m1-n1+m2-n2+m3-n3+3)./1e6, '-', 'Color', color3, 'LineWidth', 3);

errorbar(tdata(n3:m3), Tdata3(n3:m3)./1e6, st3(n3:m3)./1e6, 'o', 'Color', color3, 'MarkerSize', 8, 'MarkerFaceColor', color3, 'LineStyle', '--', 'CapSize', 3, 'LineWidth', 2);


% Adding annotations and a vertical line for treatment indication
xline(12,'--', 'm', 'LineWidth', 2.5, 'Label', '8Gy Treatment','fontsize',17, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'middle');

% Configuring axes labels and limits
xlim([7,31])
ylim([-0.11e8./1e6, 1.45e9./1e6]); % Adjust limits based on data
set(gca, 'FontSize', 19, 'LineWidth', 1.3);
grid on;
legend( 'Primary Control','Prediction','Primary 8Gy','Abscopal Control','Prediction' ,'Abscopal 8Gy ','Location', 'westoutside', 'FontSize', 20);

xlabel('\fontsize{20}Time (days)',FontWeight='bold');
ylabel('\fontsize{20}Tumor volume (mm^3)',FontWeight='bold');
%annotation('textbox', [0.05, 0.85, 0.1, 0.1], 'String', 'C', 'EdgeColor', 'none', 'FontSize', 17, 'FontWeight', 'bold');
    

function dydt = g1(t, y, x, d)
    dydt = zeros(3,1);
    t0 = d / 6;
    if d == 0
        u = 1;
    else
        u = 2 * (2.558775384581 * t0 + exp(-2.558775384581 * t0) - 1) / (2.558775384581^2 * t0^2);
    end
    S1 = exp(-0.005147638197 * d - 0.006577733464 * d^2 * u);

    dydt(1) =  0.453 * y(1) * (1 - y(1) / 1.77e9) - 0.000000044133 * y(1) * y(2) - 0.000000907400 * y(1) * y(3);
    dydt(2) = 16338.485039659030 - 0.069125930809 * y(2) - 0.000000001991 * y(1) * y(2);
    dydt(3) = 0.000000439623 - 0.394943790649 * y(3) - 0.000000815338 * y(1) * y(3);

    if (t >=12)
        dydt(1) = dydt(1) * S1 - dydt(1) * (1 - S1) * 6.12;
    end
end

function dydt = g2(t, y, x, d)
    dydt = zeros(3,1);
    t0 = d / 6;
    if d == 0
        u = 1;
    else
        u = 2 * (2.558775384581 * t0 + exp(-2.558775384581 * t0) - 1) / (2.558775384581^2 * t0^2);
    end
    S1 = exp(-0.005147638197 * d - 0.006577733464 * d^2 * u);

    dydt(1) = 0.453 * y(1) * (1 - y(1) / 1.77e9) - 0.000000044133 * y(1) * y(2) - 0.000000907400 * y(1) * y(3);
    dydt(2) = 16338.485039659030 - 0.069125930809 * y(2) - 0.000000001991 * y(1) * y(2);
    dydt(3) = 0.000000439623 - 0.394943790649 * y(3) - 0.000000815338 * y(1) * y(3);

    if (t >12)
        dydt(1) = dydt(1) - dydt(1) * (1 - S1) *6.12; %because tumor size on application day is smaller than 100
    end
end
function cost = cost_function(x, tdata, n1, n2, n3, n4, m1, m2, m3, m4, d)
    [~, Y1_common] = ode23s(@(t, y) g1(t, y, x, d(1)), tdata(n1:m1), [x(1); x(5); x(6)]);
    [~, Y2_common] = ode23s(@(t, y) g1(t, y, x, d(2)), tdata(n2:m2), [x(2); x(5); x(6)]);
    [~, Y3_common] = ode23s(@(t, y) g2(t, y, x, d(2)), tdata(n3:m3), [x(3); x(5); x(6)]);
    [~, Y4_common] = ode23s(@(t, y) g2(t, y, x, d(1)), tdata(n4:m4), [x(4); x(5); x(6)]);
 
    
    cost = [Y1_common(:, 1) + Y1_common(:, 2) + Y1_common(:, 3);
            Y2_common(:, 1) + Y2_common(:, 2) + Y2_common(:, 3);
            Y3_common(:, 1) + Y3_common(:, 2) + Y3_common(:, 3);
            Y4_common(:, 1) + Y4_common(:, 2) + Y4_common(:, 3)];

    assert(size(cost, 1) == size(tdata(n1:m1), 1) + size(tdata(n2:m2), 1) + size(tdata(n3:m3), 1) + size(tdata(n4:m4), 1), 'Inconsistent dimensions');
end
