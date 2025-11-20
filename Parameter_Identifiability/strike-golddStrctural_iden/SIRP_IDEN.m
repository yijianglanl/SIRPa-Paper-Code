% Tumor-Immune 2-compartment model for STRIKE-GOLDD
% Example model file to check structural identifiability

clear;

% 2 states
syms T I M
x = [T; I; M];

% 1 output
h = T+I+M;

% no measured input
u = [];

% 6 unknown parameters
syms x1 x2 x3 x4
p = [x1; x2; x3; x4];

f = [ 0.453*T*(1 - T/1.7e9) - 4.5e-08 *T*I - x1*T*M;
     1.65e4 - 7.41e-2*I -  1.98e-9*T*I;
      x2 - x3*M - x4*I*M ];


% initial conditions (you can specify if needed)
ics = [];         
known_ics = [1432030,94089,51640];

% save to MAT file
save('my_model2','x','p','h','f','u','ics','known_ics');
