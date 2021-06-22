function [out] = drag_forces_20m_IBP_v2(in,mdl,rb,waves,turbulence,ADCP)
%the Inertial coordinate system is has X into the flow, Z down, and Y to complete the right hand rule
%the body fixed coordinate system is attached to the CG with x to bow, y to
%stbd, and z to bottom
%bow,
u = in(1);
v = in(2);
w = in(3);
p = in(4);
p_r = in(5);
q = in(6);
r = in(7);
X = in(8);
Y = in(9);
Z = in(10);
phi = in(11);
phi_r = in(12);
theta = in(13);
psi = in(14);
Uw = in(15); %North Water Velocity
Vw = in(16); %East Water Velocity
Ww = in(17); %Down Water Velocity
dUw = in(18); %North Water Velocity Gradient
dVw = in(19); %East Water Velocity Gradient
dWw = in(20); %Down Water Velocity Gradient
t=in(21);
t_delayed = in(22);
B3_ang_rad = in(23);
B2_ang_rad = in(24);
B1_ang_rad = in(25);


%Transforms from body to inertial coordinates
L_IB = [[cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)];...
        [cos(theta)*sin(psi) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi)];...
        [-sin(theta)         sin(phi)*cos(theta)                            cos(phi)*cos(theta)                           ]];
%Transforms from inertial to body coordinates
L_BI = L_IB';
%Transform water velocity to the body fixed coordinate system
%number of discrete elements used for strip theory
N = 11; %this must be an odd number and the number of sections are N-1
S = ones(1,N)/3;
S(2:2:end-1) = 4/3;
S(3:2:end-2) = 2/3;

%body location vector
x_vec_body = mdl.CDx_body+mdl.L_body*(-0.5:1/(N-1):0.5);
y_vec_body = mdl.CDy_body*ones(1,N);
z_vec_body = mdl.CDz_body*ones(1,N);

%variable buoyancy chamber
x_vec_VBC = mdl.CDx_VBC+mdl.L_VBC*(-0.5:1/(N-1):0.5);
y_vec_VBC = mdl.CDy_VBC*ones(1,N);
z_vec_VBC = mdl.CDz_VBC*ones(1,N);

%% Wave part
if waves.H==0 %No waves
    U_wave_waves=zeros(1,3*N+3);
    V_wave_waves=zeros(1,3*N+3);
    W_wave_waves=zeros(1,3*N+3);
else %waves
    x_vec_waves=[x_vec_body,x_vec_VBC];
    y_vec_waves=[y_vec_body,y_vec_VBC];
    z_vec_waves=[z_vec_body,z_vec_VBC];
    %Inertial location vector
    X_vec_waves = X + L_IB(1,1)*x_vec_waves + L_IB(1,2)*y_vec_waves + L_IB(1,3)*z_vec_waves;
    Y_vec_waves = Y + L_IB(2,1)*x_vec_waves + L_IB(2,2)*y_vec_waves + L_IB(2,3)*z_vec_waves;
    Z_vec_waves = Z + L_IB(3,1)*x_vec_waves + L_IB(3,2)*y_vec_waves + L_IB(3,3)*z_vec_waves;

    [w2,X_vec_waves2]=meshgrid(waves.w2,X_vec_waves);
    [k2,Y_vec_waves2]=meshgrid(waves.k2,Y_vec_waves);
    [teta2,Z_vec_waves2]=meshgrid(waves.teta2,Z_vec_waves);
    [H,~]=meshgrid(waves.H,ones(1,3*N+3));
    [pri,~]=meshgrid(waves.pri,ones(1,3*N+3));

    UV_wave_waves=H.*waves.g.*k2./w2.*cosh(k2.*(waves.h-Z_vec_waves2))./cosh(k2.*waves.h).*sin(k2.*(X_vec_waves2.*cos(teta2)+Y_vec_waves2.*sin(teta2))-w2.*t+pri); 
    U_wave_waves=sum(cos(teta2).*UV_wave_waves,2)';%calculates the forward velocity of the particle of water
    V_wave_waves=sum(sin(teta2).*UV_wave_waves,2)';%calculates the forward velocity of the particle of water
    W_wave_waves=sum(H.*waves.g.*k2./w2.*sinh(k2.*(waves.h-Z_vec_waves2))./cosh(k2.*waves.h).*cos(k2.*(X_vec_waves2.*cos(teta2)+Y_vec_waves2.*sin(teta2))-w2.*t+pri),2)'; %calculates the downward velocity of the particle of water
end

U_wave_body=U_wave_waves(1:N); %North Component of Wave Induced Velocity 
U_wave_VBC=U_wave_waves(N+1:2*N); 

V_wave_body=V_wave_waves(1:N); %East Component of Wave Induced Velocity
V_wave_VBC=V_wave_waves(N+1:2*N);

W_wave_body=W_wave_waves(1:N); %Down Component of Wave Induced Velocity 
W_wave_VBC=W_wave_waves(N+1:2*N);

%% Calculate the drag on the body and Variable Buoyancy Chamber (VBC)
%****use strip theory for a cross flow
%Inertial location vector
Z_vec_body = Z + L_IB(3,1)*x_vec_body + L_IB(3,2)*y_vec_body + L_IB(3,3)*z_vec_body;
Z_vec_VBC = Z + L_IB(3,1)*x_vec_VBC + L_IB(3,2)*y_vec_VBC + L_IB(3,3)*z_vec_VBC;

%velocity vector
u_vec_body = u+q*z_vec_body-r*y_vec_body;
v_vec_body = v+r*x_vec_body-p*z_vec_body;
w_vec_body = w+p*y_vec_body-q*x_vec_body;
u_vec_VBC = u+q*z_vec_VBC-r*y_vec_VBC;
v_vec_VBC = v+r*x_vec_VBC-p*z_vec_VBC;
w_vec_VBC = w+p*y_vec_VBC-q*x_vec_VBC;

% calculate the turbulence on the center of the body and VBC
u_1_body=sum(turbulence.mag_x(:,1).*sin(2*pi*(turbulence.fm'*ones(1,1))*t+turbulence.ang_x(:,1))); %1st column of turbulence velocity in frequency domain used
v_body=sum(turbulence.mag_y(:,1).*sin(2*pi*(turbulence.fm'*ones(1,1))*t+turbulence.ang_y(:,1)));
w_body=sum(turbulence.mag_z(:,1).*sin(2*pi*(turbulence.fm'*ones(1,1))*t+turbulence.ang_z(:,1)));
u_body=u_1_body+turbulence.r_uv*v_body+turbulence.r_uw*w_body;

u_1_VBC=sum(turbulence.mag_x(:,end).*sin(2*pi*(turbulence.fm'*ones(1,1))*t+turbulence.ang_x(:,end))); %Last column of turbulence velocity in frequency domain used
v_VBC=sum(turbulence.mag_y(:,end).*sin(2*pi*(turbulence.fm'*ones(1,1))*t+turbulence.ang_y(:,end)));
w_VBC=sum(turbulence.mag_z(:,end).*sin(2*pi*(turbulence.fm'*ones(1,1))*t+turbulence.ang_z(:,end)));
u_VBC=u_1_VBC+turbulence.r_uv*v_VBC+turbulence.r_uw*w_VBC;

ADCP.inflowSelection = 1;
%water current velocity vector
%InflowSelection=1 if using an imported dataset else use uniform values.
if ADCP.inflowSelection == 1 
    time_in = ADCP.tStart+datenum([0,0,0,0,0,t]);
    Uw_interp_body_vec = interp2(ADCP.z,ADCP.t,ADCP.u',Z_vec_body(:),time_in);
    Vw_interp_body_vec = interp2(ADCP.z,ADCP.t,ADCP.v',Z_vec_body(:),time_in);
    Ww_interp_body_vec = interp2(ADCP.z,ADCP.t,ADCP.w',Z_vec_body(:),time_in);
    Uw_interp_body = reshape(Uw_interp_body_vec,size(Z_vec_body));
    Vw_interp_body = reshape(Vw_interp_body_vec,size(Z_vec_body));
    Ww_interp_body = reshape(Ww_interp_body_vec,size(Z_vec_body));
    Uw_vec_body = Uw_interp_body+u_body;
    Vw_vec_body = Vw_interp_body+v_body;
    Ww_vec_body = Ww_interp_body+w_body;
    Uw_interp_VBC_vec = interp2(ADCP.z,ADCP.t,ADCP.u',Z_vec_VBC(:),time_in);
    Vw_interp_VBC_vec = interp2(ADCP.z,ADCP.t,ADCP.v',Z_vec_VBC(:),time_in);
    Ww_interp_VBC_vec = interp2(ADCP.z,ADCP.t,ADCP.w',Z_vec_VBC(:),time_in);
    Uw_interp_VBC = reshape(Uw_interp_VBC_vec,size(Z_vec_VBC));
    Vw_interp_VBC = reshape(Vw_interp_VBC_vec,size(Z_vec_VBC));
    Ww_interp_VBC = reshape(Ww_interp_VBC_vec,size(Z_vec_VBC));
    Uw_vec_VBC = Uw_interp_VBC+u_VBC;
    Vw_vec_VBC = Vw_interp_VBC+v_VBC;
    Ww_vec_VBC = Ww_interp_VBC+w_VBC;
else
    Uw_vec_body = Uw+Z_vec_body*dUw+u_body; %North current + turbulence - replacing first 2 terms
    Vw_vec_body = Vw+Z_vec_body*dVw+v_body; %East current + turbulence
    Ww_vec_body = Ww+Z_vec_body*dWw+w_body; %Down current + turbulence
    Uw_vec_VBC = Uw+Z_vec_VBC*dUw+u_VBC;%turbulence taken into account
    Vw_vec_VBC = Vw+Z_vec_VBC*dVw+v_VBC;
    Ww_vec_VBC = Ww+Z_vec_VBC*dWw+w_VBC;
end

uw_vec_body = L_BI(1,1)*Uw_vec_body + L_BI(1,2)*Vw_vec_body + L_BI(1,3)*Ww_vec_body; % body current + turbulence
vw_vec_body = L_BI(2,1)*Uw_vec_body + L_BI(2,2)*Vw_vec_body + L_BI(2,3)*Ww_vec_body;
ww_vec_body = L_BI(3,1)*Uw_vec_body + L_BI(3,2)*Vw_vec_body + L_BI(3,3)*Ww_vec_body;
uw_vec_VBC = L_BI(1,1)*Uw_vec_VBC + L_BI(1,2)*Vw_vec_VBC + L_BI(1,3)*Ww_vec_VBC;
vw_vec_VBC = L_BI(2,1)*Uw_vec_VBC + L_BI(2,2)*Vw_vec_VBC + L_BI(2,3)*Ww_vec_VBC;
ww_vec_VBC = L_BI(3,1)*Uw_vec_VBC + L_BI(3,2)*Vw_vec_VBC + L_BI(3,3)*Ww_vec_VBC;

%water wave velocity vector for body and VBC
u_wave_body = L_BI(1,1)*U_wave_body + L_BI(1,2)*V_wave_body + L_BI(1,3)*W_wave_body;
v_wave_body = L_BI(2,1)*U_wave_body + L_BI(2,2)*V_wave_body + L_BI(2,3)*W_wave_body;
w_wave_body = L_BI(3,1)*U_wave_body + L_BI(3,2)*V_wave_body + L_BI(3,3)*W_wave_body;
u_wave_VBC = L_BI(1,1)*U_wave_VBC + L_BI(1,2)*V_wave_VBC + L_BI(1,3)*W_wave_VBC;
v_wave_VBC = L_BI(2,1)*U_wave_VBC + L_BI(2,2)*V_wave_VBC + L_BI(2,3)*W_wave_VBC;
w_wave_VBC = L_BI(3,1)*U_wave_VBC + L_BI(3,2)*V_wave_VBC + L_BI(3,3)*W_wave_VBC;

%relative velocity vector for body and VBC
u_rel_body = u_vec_body-uw_vec_body-u_wave_body; %body-current-waves-turbulence
v_rel_body = v_vec_body-vw_vec_body-v_wave_body;
w_rel_body = w_vec_body-ww_vec_body-w_wave_body;
u_rel_VBC = u_vec_VBC-uw_vec_VBC-u_wave_VBC;
v_rel_VBC = v_vec_VBC-vw_vec_VBC-v_wave_VBC;
w_rel_VBC = w_vec_VBC-ww_vec_VBC-w_wave_VBC;

%calculate the forces on the body and VBC
fx_body_vec = -0.5*mdl.rho*mdl.Cdx_body*mdl.D_body^2*pi/4*u_rel_body.*sqrt(u_rel_body.^2+v_rel_body.^2+w_rel_body.^2);%
fx_body = sum(fx_body_vec.*S);
fy_body_vec = -0.5*mdl.rho*mdl.Cdy_body*mdl.D_body*mdl.L_body/N*v_rel_body.*sqrt(u_rel_body.^2+v_rel_body.^2+w_rel_body.^2);%
fy_body = sum(fy_body_vec.*S);
fz_body_vec = -0.5*mdl.rho*mdl.Cdz_body*mdl.D_body*mdl.L_body/N*w_rel_body.*sqrt(u_rel_body.^2+v_rel_body.^2+w_rel_body.^2);%
fz_body = sum(fz_body_vec.*S);

fx_VBC_vec = -0.5*mdl.rho*mdl.Cdx_VBC*mdl.D_VBC^2*pi/4/N*u_rel_VBC.*sqrt(u_rel_VBC.^2+v_rel_VBC.^2+w_rel_VBC.^2);%
fx_VBC = sum(fx_VBC_vec.*S);
fy_VBC_vec = -0.5*mdl.rho*mdl.Cdy_VBC*mdl.D_VBC*mdl.L_VBC/N*v_rel_VBC.*sqrt(u_rel_VBC.^2+v_rel_VBC.^2+w_rel_VBC.^2);%
fy_VBC = sum(fy_VBC_vec.*S);
fz_VBC_vec = -0.5*mdl.rho*mdl.Cdz_VBC*mdl.D_VBC*mdl.L_VBC/N*w_rel_VBC.*sqrt(u_rel_VBC.^2+v_rel_VBC.^2+w_rel_VBC.^2);%
fz_VBC = sum(fz_VBC_vec.*S);

%calculate the moments on the turbine body and VBC
Mx_body = sum((fz_body_vec.*y_vec_body-fy_body_vec.*z_vec_body).*S);
My_body = sum((fx_body_vec.*z_vec_body-fz_body_vec.*x_vec_body).*S);
Mz_body = sum((fy_body_vec.*x_vec_body-fx_body_vec.*y_vec_body).*S);
Mx_VBC = sum((fz_VBC_vec.*y_vec_VBC-fy_VBC_vec.*z_vec_VBC).*S);
My_VBC = sum((fx_VBC_vec.*z_vec_VBC-fz_VBC_vec.*x_vec_VBC).*S);
Mz_VBC = sum((fy_VBC_vec.*x_vec_VBC-fx_VBC_vec.*y_vec_VBC).*S);

%% Calculate the rotor forces
%in(1:12)

[fx_rotor, fy_rotor, fz_rotor, Mx_rotor, My_rotor, Mz_rotor, Mx_s, My_s, ...
    Mz_s, MA_blade1, MA_blade2, MA_blade3, fx_oneblade] = ...
    rotor_forces_20mD_turbulence_v3([u,v,w,p,q,r,X,Y,Z, ...
    phi,theta,psi,t,t_delayed],[Uw Vw Ww dUw dVw dWw],[p_r*60/2/pi phi_r],...
    mdl,rb,waves,turbulence,ADCP,[B3_ang_rad B2_ang_rad B1_ang_rad]);

%Here are the total forces and moments
fx = fx_body+fx_VBC+fx_rotor;
fy = fy_body+fy_VBC+fy_rotor;
fz = fz_body+fz_VBC+fz_rotor;
Mx_b = Mx_body+Mx_VBC;
Mx_r = Mx_rotor;
My = My_body+My_VBC+My_rotor;
Mz = Mz_body+Mz_VBC+Mz_rotor;

%[fx fy fz Mx_b Mx_r My Mz fx_rotor fy_rotor fz_rotor Mx_s My_s Mz_s fx_mast fy_mast fz_mast Mx_mast My_mast Mz_mast MA MT]
%[U_wave_aquadopp V_wave_aquadopp W_wave_aquadopp]
out = [fx fy fz Mx_b Mx_r My Mz fx_rotor fy_rotor fz_rotor MA_blade3 MA_blade2 MA_blade1 fx_oneblade];