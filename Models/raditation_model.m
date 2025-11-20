function dydt = raditation_model(t, y, x, d)
   dydt = zeros(2,1);
   x(3)=2.56;
   t0 = d/1.2;
   if (d == 0)
       u = 1;
   else
       u = 2 * (x(3)*t0 + exp(-x(3)*t0) - 1) / ((x(3)^2)*t0^2);
   end
      
   S = exp(-5.15e-3*d - 6.58e-3*(d^2)*u);
 
   dydt(1) = 0.453 *y(1)*(1 - y(1)/1.77e9) -4.5e-08 *y(1)*y(2);
   dydt(2) = 1.65e4-7.41e-2*y(2)- 1.98e-9*y(1)*y(2);
    
   if (t >= d && d > 0)
       dydt(1) = dydt(1) * S;
   end
end
