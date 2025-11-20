% Tumor-Immune 2-compartment model for STRIKE-GOLDD
clear;

% States (column vector)
syms T I real
x = [T; I];

% Output
h = T + I;

% No measured inputs
u = [];

% Parameters (9 unknowns)
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 real
p = [x1; x2; x3; x4; x5; x6; x7; x8; x9];

% --- Your dose-specific auxiliary (rename to avoid clashing with 'u') ---
d = sym(8);                 % if you want dose fixed at 8
t0 = d/1.2;
ud = 2*( x3*t0 + exp(-x3*t0) - 1 ) / ( (x3^2) * t0^2 );   % renamed from 'u'

% Treatment scaling S1(d)
S1 = exp(-x1*d - x2*d^2*ud);

% Dynamics
f = [ S1*( x4*T*(1 - T/x5) - x6*T*I );
      x7 - x9*I - x8*T*I ];

% Initial conditions
ics = [];                   % unknown ICs (none specified)
known_ics = [1e5; 1e4];     % make it a column to match x

% Save model
save('my_modelWT','x','p','h','f','u','ics','known_ics');
