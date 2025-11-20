%MC38_SIPRA_KO

function dydt = MC38_SIRPA_KO(~,y,x)
    dydt = zeros(3,1);
    dydt(1) = 0.453 *y(1)*(1 - y(1)/1.77e9) -4.5e-08 *y(1)*y(2)-x(1)*y(1)*y(3);
    dydt(2) = 1.65e4-7.41e-2*y(2)- 1.98e-9*y(1)*y(2);
    dydt(3) = x(2) - x(3)*y(3) - x(4)*y(1)*y(3);
end
