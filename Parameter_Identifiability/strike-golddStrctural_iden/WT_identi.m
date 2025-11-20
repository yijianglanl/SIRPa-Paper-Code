% Tumor-Immune 2-compartment model for STRIKE-GOLDD
% Example model file to check structural identifiability

clear;

% 2 states
syms T I
x = [T; I];

% 1 output
h = T+I;

% no measured input
u = [];

% 6 unknown parameters
syms x1 x2 x3 x4 x5 x6
p = [x1; x2; x3; x4; x5; x6];

% dynamic equations (copying your ODE definitions)
f = [ x1*T*(1-(T/x2)) - x3*T*I ;
      x4 - x5*I - x6*T*I ];

% initial conditions (you can specify if needed)
ics = [];         
known_ics = [1e5,1e4];

% save to MAT file
save('my_model1','x','p','h','f','u','ics','known_ics');
