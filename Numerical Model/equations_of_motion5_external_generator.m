function [out] = equations_of_motion5_external_generator(in,mdl)

u = in(1); %surge velocity of the turbine
v = in(2); %sway velocity of the turbine
w = in(3); %heave velocity of the turbine
p_b = in(4); %rotational velocity of the body about x 
p_r = in(5); %rotational velocity of the rotor about x 
q = in(6); %rotational velocity of the turbine about y 
r = in(7); %rotational velocity of ther turbine about z
%X = in(8); %not used
%Y = in(9);    %not used
%Z = in(10); %not used
phi = in(11); %roll of the body
%phi_r = in(12); %roll of the rotor
theta = in(13); %pitch of the turbine
psi = in(14);   %yaw of the turbine
fx = in(15); %hydrodynamic and hydrostatic forces on the turbine in x direction
fy = in(16); %hydrodynamic and hydrostatic forces on the turbine in y direction
fz = in(17); %hydrodynamic and hydrostatic forces on the turbine in z direction
Mx_b = in(18); %net hydrodynamic, hydrostatic, and internal moment on the main body (not the entire turbine) about the x axis
Mx_r = in(19); %net hydrodynamic, hydrostatic, and internal moment on the rotor about the x axis
My = in(20); %net hydrodynamic and hydrostatic moment on the turbine about the y axis
Mz = in(21); %net hydrodynamic and hydrostatic moment on the turbine about the z axis

FX_c = in(22);
FY_c = in(23);
FZ_c = in(24);

Mx_s = in(25); %in(25) is the electromechanical torque placed on the rotor about the x-axis -> an equal but opisite moment should be placed on the main body

%coordinate transformations
L_IB = [[cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)];...
        [cos(theta)*sin(psi) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi)];...
        [-sin(theta)         sin(phi)*cos(theta)                            cos(phi)*cos(theta)                           ]]; 
%Transforms from inertial to body coordinates
L_BI = L_IB';
%%%%
T_IB = [[1, sin(phi)*tan(theta), cos(phi)*tan(theta)];[0, cos(phi), -sin(phi)];[0, sin(phi)/cos(theta), cos(phi)/cos(theta)]];

%calculate cable forces in the body fixed frame
f_c = L_BI*[FX_c FY_c FZ_c]';
fx_t = fx+f_c(1);
fy_t = fy+f_c(2);
fz_t = fz+f_c(3);

%calculate the moments
M_c = cross([mdl.CAPx mdl.CAPy mdl.CAPz],f_c');

Mx_b_t = Mx_b + M_c(1); %these are only the moments on the body and not the entire turbine
Mx_r_t = Mx_r;
My_t = My + M_c(2);
Mz_t = Mz + M_c(3);

%**************************
%this section has been updated since Nicolas
%**************************
m_v = 2*mdl.mass_T; %vertual mass of the entire turbine
m_b = 2*mdl.mass_b; %vertual mass of the body
m_r = 2*mdl.mass_r; %vertual mass of the rotor section

% Ix = 2*mdl.Ixx; %vertual mass moment of inertia of the turbine about x 
% Iy = 2*mdl.Iyy; %vertual mass moment of inertia of the turbine about y 
% Iz = 2*mdl.Izz; %vertual mass moment of inertia of the turbine about z 
% Ixz = 2*mdl.Ixz; %vertual mass product of inertia of the turbine about x z 

Ix_b = 2*mdl.Ixx_b; %vertual mass moment of inertia of the body about x 
Iy_b = 2*mdl.Iyy_b; %vertual mass moment of inertia of the body about y 
Iz_b = 2*mdl.Izz_b; %vertual mass moment of inertia of the body about z 
Ixz_b = 2*mdl.Ixz_b; %vertual mass product of inertia of the body about x z 

Ix_r = 2*mdl.Ixx_r; %vertual mass moment of inertia of the rotor section about x 
Iy_r = 2*mdl.Iyy_r; %vertual mass moment of inertia of the rotor section about y 
Iz_r = 2*mdl.Izz_r; %vertual mass moment of inertia of the rotor section about z
%Ixz_r = 2*mdl.Ixz_r; %vertual mass product of inertia of the rotor about x z 
mat_inverse = mdl.mat_inv;

x_s = mdl.CMx; %center of mass of the entire system in x
%y_s = mdl.CMy; %center of mass of the entire system in y
% z_s = mdl.CMz; %center of mass of the entire system in z
x_b = mdl.CMx_b;
%y_b = mdl.CMy_b;
z_b = mdl.CMz_b;
x_r = mdl.CMx_r;
%y_r = mdl.CMy_r;
% z_r = mdl.CMz_r;

%Calculate the shaft moment for a fixed RPM

% r_dot = 0;
% for ct=1:7
%     A_var = q*r*(Iz_b-Iy_b)-r_dot*Ixz_b-p_b*q*Ixz_b+w*m_b*z_b*p_b;
%     B_var = q*r*(Iz_r-Iy_r)-r_dot*Ixz_r-p_r*q*Ixz_r+w*m_b*z_r*p_r;
%     Mx_s = (Mx_r_t*Ix_b-Mx_b_t*Ix_r+A_var*Ix_r-B_var*Ix_b)/(Ix_b+Ix_r);

%uncoupled equation of motion
%p_r_dot = (Mx_r_t - Mx_s - q*r*(Iz_r-Iy_r));%/Ix_r; %this is no longer the moment of inertia but instead it is what is fed into the generator model!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
p_r_dot_ex_gen = (Mx_r_t - q*r*(Iz_r-Iy_r));%/Ix_r;

acc6 = mat_inverse*...
[fx_t + m_v*(v*r-w*q) + m_v*x_s*(q^2+r^2)-m_b*z_b*p_b*r;...
 fy_t - m_v*u*r + w*(m_b*p_b+m_r*p_r) - m_b*z_b*q*r - m_b*x_b*q*p_b - m_r*x_r*q*p_r;...
 fz_t + m_v*u*q - v*(m_b*p_b+m_r*p_r) + m_b*z_b*(p_b^2+q^2) - m_b*x_b*r*p_b - m_r*x_r*r*p_r;...
 Mx_b_t - Mx_s - q*r*(Iz_b-Iy_b) + Ixz_b*p_b*q - m_b*z_b*(w*p_b-u*r);... % I changed Mx_b_t + Mx_s to Mx_b_t - Mx_s on 11/19/2015
 My_t - r*p_b*(Ix_b-Iz_b) - r*p_r*(Ix_r-Iz_r) - Ixz_b*(p_b^2-r^2) + m_b*z_b*(v*r-w*q) - m_v*x_s*u*q + m_b*x_b*v*p_b + m_r*x_r*v*p_r;...
 Mz_t - q*p_b*(Iy_b-Ix_b) - q*p_r*(Iy_r-Ix_r) - Ixz_b*r*q - m_v*x_s*u*r + m_b*x_b*w*p_b + m_r*x_r*w*p_r];

u_dot = acc6(1);
v_dot = acc6(2);
w_dot = acc6(3);
p_b_dot = acc6(4);
q_dot = acc6(5);
r_dot = acc6(6);
    
%Linear velocities
VEL = L_IB*[u; v; w];
X_dot = VEL(1);
Y_dot = VEL(2);
Z_dot = VEL(3);
%Euler velocities
ANG_VEL_b = T_IB*[p_b; q; r];
phi_b_dot = ANG_VEL_b(1);
theta_dot = ANG_VEL_b(2);
psi_dot = ANG_VEL_b(3);

%****************************
%end updated section
%****************************



out = [u_dot v_dot w_dot p_b_dot q_dot r_dot X_dot Y_dot Z_dot phi_b_dot theta_dot psi_dot p_r_dot_ex_gen];