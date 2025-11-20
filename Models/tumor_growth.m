%tumor growth Model
function dydt = tumor_growth(~,y,x)

dydt = zeros(2,1);
dydt(1) = x(1)*y(1)*(1 - y(1)/x(2)) - x(3)*y(1)*y(2);
dydt(2) = x(4) - x(5)*y(2) - x(6)*y(1)*y(2);

end