% Tumor-Immune 2-compartment model for STRIKE-GOLDD
clear;

% States (column vector)
syms T I M real
x = [T; I;M];

% Output
h = T + I + M;

% No measured inputs
u = [];

% Parameters (9 unknowns)
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15  real
p = [x1; x2; x3; x4; x5; x6; x7; x8; x9;x10; x11; x12; x13; x14; x15];

% --- Your dose-specific auxiliary (rename to avoid clashing with 'u') ---
d = sym(8);                 % if you want dose fixed at 8
t0 = d/1.2;
ud = 2*( x3*t0 + exp(-x3*t0) - 1 ) / ( (x3^2) * t0^2 );   % renamed from 'u'

% Treatment scaling S1(d)
S1 = exp(-x1*d - x2*d^2*ud);
%16338.485039659030 -   0.069125930809* I - 0.000000001991* I * T
%0.000000439623   -0.394943790649 *y(3)-0.000000815338*y(1)*y(3)
% Dynamics
f = [ S1*( x4*T*(1 - T/x5) - x6*T*I-x10*T*M )-(1-S1)*x15*( x4*T*(1 - T/x5) - x6*T*I-x10*T*M );
      x7 - x9*I - x8*T*I;
      x11 - x12*M - x13*T*M];

% Initial conditions
ics = [];                   % unknown ICs (none specified)
known_ics = [1e5; 1e4; 1e4];     % make it a column to match x

% Save model
save('my_modelWT','x','p','h','f','u','ics','known_ics');
