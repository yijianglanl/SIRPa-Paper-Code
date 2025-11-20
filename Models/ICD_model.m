%ICD_model
function dydt = ICD_model(t, y, x,d)
    dydt = zeros(3,1);
        t0 = d/1.2;
        if (d==0)
            u=1;
        else
            u = 2 * (2.56 *t0 + exp(-2.56 *t0) - 1) / ((  2.56^2)*t0^2);
        end
       S=exp(-5.15e-3 *d -6.58e-3*d*d*u);
    if d == 4
        A = x(1);
    elseif d == 8
         A = x(2);
    elseif d == 15
         A = x(3);
    else
        A=0;  
    end
    dydt(1) = 0.453 *y(1)*(1 - y(1)/1.77e9) -4.5e-08 *y(1)*y(2)-9.07e-7*y(1)*y(3);
    dydt(2) = 1.65e4-7.41e-2*y(2)- 1.98e-9*y(1)*y(2);
    dydt(3) =4.39e-7   -3.949e-1 *y(3)-8.15338e-7*y(1)*y(3);
    
    if (t >= 8)
        dydt(1) = dydt(1)*S-A*(1-S)*dydt(1); 
    else
        dydt(1) = dydt(1); 
    end
    end