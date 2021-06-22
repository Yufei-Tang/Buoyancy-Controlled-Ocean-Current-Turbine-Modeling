function [out] = attachment_point(in,mdl)
u = in(1);
v = in(2);
w = in(3);
p = in(4);
%p_r = in(5)
q = in(6);
r = in(7);
X = in(8);
Y = in(9);
Z = in(10);
phi = in(11);
%phi_r = in(12);
theta = in(13);
psi = in(14);

%Transforms from body to inertial coordinates
L_IB = [[cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)];...
        [cos(theta)*sin(psi) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi)];...
        [-sin(theta)         sin(phi)*cos(theta)                            cos(phi)*cos(theta)                           ]];

%Calculate the position of the nose in the inertial frame
POS = [X Y Z]'+ L_IB*[mdl.CAPx mdl.CAPy mdl.CAPz]';
vel = [u v w]+ cross([p q r],[mdl.CAPx mdl.CAPy mdl.CAPz]);
VEL = L_IB*vel';
 
out = [VEL(1) POS(1) VEL(2) POS(2) VEL(3) POS(3)];