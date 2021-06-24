function [torque_EM] = controller_v2(RPM,C_pmax,TSR_max)
% place all of your equations here that take you from your inputs to the 
% electromechanical torque on the rotor shaft that apposses rotor rotation
%According to the C_p calculation code my C_pmax and TSR that yields this
%value is the following
%C_pmax = 0.4618;
%TSR_max = 3.9516;
mdl.rho = 1025.2;
R = 10;
% Calculate cross-sectional area of rotor
A = R*R*pi;
K = 0.5*mdl.rho*A*R*R*R*C_pmax/(TSR_max^3);
% Convert RPM to rad/s
w = RPM*2*pi/60;
torque_EM = K*w*w;