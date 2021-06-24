function [fx_t, fy_t, fz_t, Mx, My, Mz, Mx_s, My_s, Mz_s, MA_blade1, MA_blade2, MA_blade3, fx_oneblade] = rotor_forces_20mD_turbulence_v3(states,U_vec,rb_states,mdl,rb,waves,turbulence,ADCP,root_twistB3B2B1)
%this version changes the coordinate system so that the attack and twist
%angles are measured with respect to the rottor axis
N_mesh = rb.N_mesh;
B3_ang_rad = root_twistB3B2B1(1);
B2_ang_rad = root_twistB3B2B1(2); 
B1_ang_rad = root_twistB3B2B1(3);
root_twistB1B2B3 = [B1_ang_rad B2_ang_rad B3_ang_rad];

%here are the global variables that should be passed
global W_n_b W_t_b W_n_m W_t_m W_0_n_m W_0_t_m W_int_n_m W_int_t_m W_qs_n_m W_qs_t_m time_delayed_old root_twist_m_global ang_global;
global W_n_b_old W_t_b_old W_n_m_old W_t_m_old W_0_n_m_old W_0_t_m_old W_int_n_m_old W_int_t_m_old W_qs_n_m_old W_qs_t_m_old root_twist_m_global_old ang_global_old;

%W_n = normal wake on the actual rotor

%U_vec = Uw Vw Ww dUw/dZ dVw/dZ dWw/dZ

%these rotor blade states should be lumped into the system states once the
%simulation is updated for the new equations of motion
RPM = rb_states(1);
ang = rb_states(2); %this as also theta for blade 1

%Turbine system states for 6DOF EOM 
u = states(1);
v = states(2);
w = states(3);
p = states(4);
q = states(5);
r = states(6);
X = states(7);
Y = states(8);
Z = states(9);
phi = states(10);
theta = states(11);
psi = states(12);
t = states(13);         %this is the actual time for the system (input from clock)
time_delayed=states(14);%this is the time of the last actual time step (clock with time delay)

%time since the last time step
delta_t = t - time_delayed;

if time_delayed_old==time_delayed %This happens when an actual time step was not taken
    W_n_b=W_n_b_old;
    W_t_b=W_t_b_old;
    W_n_m=W_n_m_old;
    W_t_m=W_t_m_old;    
    W_0_n_m=W_0_n_m_old;
    W_0_t_m=W_0_t_m_old;
    W_int_n_m=W_int_n_m_old;
    W_int_t_m=W_int_t_m_old;
    W_qs_n_m=W_qs_n_m_old;
    W_qs_t_m=W_qs_t_m_old;
    root_twist_m_global = root_twist_m_global_old;
else %this is true when an actual time step is taken
    W_n_b_old=W_n_b;
    W_t_b_old=W_t_b;
    W_n_m_old=W_n_m;
    W_t_m_old=W_t_m;
    W_0_n_m_old=W_0_n_m;
    W_0_t_m_old=W_0_t_m;
    W_int_n_m_old=W_int_n_m;
    W_int_t_m_old=W_int_t_m;
    W_qs_n_m_old=W_qs_n_m;
    W_qs_t_m_old=W_qs_t_m;
    time_delayed_old=time_delayed;
    root_twist_m_global_old = root_twist_m_global;
    ang_global_old = ang_global;
end

%Transforms from body to inertial coordinates
L_IB = [[cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)];...
    [cos(theta)*sin(psi) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi)];...
    [-sin(theta)         sin(phi)*cos(theta)                            cos(phi)*cos(theta)                           ]];
%Transforms from inertial to body coordinates
L_BI = L_IB';

%*******************************CONSTANTS*******************************
%These should be inputs (same for all 3 blades)
r_sec = rb.radius; %radius for the different sections
c_sec = rb.chord; %cord vector for the different sections
b_sec = rb.pre_twist*pi/180; %angle of rotor blade with respect to rotor plane

%these create matricies that have the same values for each row as the original vectors and has the same number of rows as the number of rotor blades  
%matricies on the rotor blades
c_b = meshgrid(c_sec,zeros(1,rb.N_blades))'; %cord matrix for rotor blades
r_b = meshgrid(r_sec,zeros(1,rb.N_blades))'; %radius matrix for rotor blades
b_b = meshgrid(b_sec,zeros(1,rb.N_blades))'+ones(1,25)'*[B1_ang_rad B2_ang_rad B3_ang_rad]; %angle wrt rotor plane matrix for rotor blades (4pi/3+ang 2pi/3+ang ang)
% this updates the root twist global variable for the three rotor blades
ang_global = ang+[4*pi/3 2*pi/3 0];%2*pi*(1/rb.N_blades:1/rb.N_blades:1);
% this updates the momentum mesh twist angles based on the blade angles
% that passed the mesh during the last time step
%this is the index of the mesh angles that were passed and the coresponding
%blade that passed it
ang_change = ang_global-ang_global_old;
ang_old = atan2(sin(ang_global_old),cos(ang_global_old));
ang_new = ang_old + ang_change;
ang_mesh = [2*pi*(1/N_mesh:1/N_mesh:1)-2*pi 2*pi*(1/N_mesh:1/N_mesh:1) 2*pi+2*pi*(1/N_mesh:1/N_mesh:1)];
I1_all = find(ang_mesh < ang_new(1) & ang_mesh > ang_old(1));
I1 = I1_all - floor((I1_all-0.000001)/N_mesh)*N_mesh;
I2_all = find(ang_mesh < ang_new(2) & ang_mesh > ang_old(2));
I2 = I2_all - floor((I2_all-0.000001)/N_mesh)*N_mesh;
I3_all = find(ang_mesh < ang_new(3) & ang_mesh > ang_old(3));
I3 = I3_all - floor((I3_all-0.000001)/N_mesh)*N_mesh;
index_passed_m = [I1 I2 I3];
index_passed_b = [I1*0+1 I2*0+2 I3*0+3];
%this updates the mesh angles that were passed
if length(index_passed_m) > 0.5
    root_twist_m_global(index_passed_m) = root_twistB1B2B3(index_passed_b);
end
%this updates the 
%matricies on the blade momentum mesh
c_m = meshgrid(c_sec,zeros(1,N_mesh))'; %cord matrix for momentum model
r_m = meshgrid(r_sec,zeros(1,N_mesh))'; %radius matrix for momentum model
b_m = meshgrid(b_sec,zeros(1,N_mesh))'+ones(1,25)'*root_twist_m_global; %angle wrt rotor plane matrix for momentum model %%%%%%%UPDATE !!!!!!!!!!!
%S_m = meshgrid(rb.S,zeros(1,N_mesh))'; %simsons vector for momentum model
%***************************************************************
%calculate the location of the blade elements with respect to the rotor hub
y_rh = r_sec*sin(ang_global); %y location of rotor blade wrt rotor hub
z_rh = -r_sec*cos(ang_global);%z location of rotor blade wrt rotor hub

%calculate the location of the mesh grid with respect to the rotor hub
%******************CONSTANT VALUES - MOVE TO THE CONSTANTS SECTION*********
y_rh_m = r_sec*sin(2*pi*(1/N_mesh:1/N_mesh:1));  %y location of the meshgrid with respect to the rotor hub
z_rh_m = -r_sec*cos(2*pi*(1/N_mesh:1/N_mesh:1)); %z location of the meshgrid with respect to the rotor hub
%**************************************************************************

% body fixed location of the sections of the blades and momentum mesh
%sections of the blade
x_b = mdl.CDx_rotor+0*y_rh; %x-location of the sections of the blades
y_b = mdl.CDy_rotor+y_rh; %y-location of the sections of the blades
z_b = mdl.CDz_rotor+z_rh; %z-location of the sections of the blades

%sections of the mesh
%******(CONSTANT VALUES - MOVE TO THE CONSTANTS SECTION)***************
x_m = mdl.CDx_rotor+0*y_rh_m; %x-location of the sections of the momentum mesh
y_m = mdl.CDy_rotor+y_rh_m;   %y-location of the sections of the momentum mesh
z_m = mdl.CDz_rotor+z_rh_m;   %z-location of the sections of the momentum mesh
%**********************************************************************

%calculate the inertial location of the section of the blades
X_b = X + L_IB(1,1)*x_b + L_IB(1,2)*y_b + L_IB(1,3)*z_b;
Y_b = Y + L_IB(2,1)*x_b + L_IB(2,2)*y_b + L_IB(2,3)*z_b;
Z_b = Z + L_IB(3,1)*x_b + L_IB(3,2)*y_b + L_IB(3,3)*z_b;

%calculate the inertial location of the section of the momentum mesh
X_m = X + L_IB(1,1)*x_m + L_IB(1,2)*y_m + L_IB(1,3)*z_m;
Y_m = Y + L_IB(2,1)*x_m + L_IB(2,2)*y_m + L_IB(2,3)*z_m;
Z_m = Z + L_IB(3,1)*x_m + L_IB(3,2)*y_m + L_IB(3,3)*z_m;

%*****************************************************************
%This is the velocities of the blade sections V_rot
%*****************************************************************
%matrix of rotor blade angles
ang_b = ones(length(Z_b),1)*ang_global;
ang_m = ones(length(Z_b),1)*2*pi*(1/N_mesh:1/N_mesh:1);

% calculate the velocity of the blade sections (9) V_rot
%       CG : vel from turb rot : vel from blade rot
u_b = u  +  q*z_b-r*y_b  +  0;
v_b = v  +  r*x_b-p*z_b  -  RPM/60*2*pi*z_rh;
w_b = w  +  p*y_b-q*x_b  +  RPM/60*2*pi*y_rh;
n_b = -u_b; %affect of rotor velocity on normal relative water velocity in the x direction
t_b = -w_b.*sin(ang_b) - v_b.*cos(ang_b); %affect of rotor velocity
%on tangental relative water velocity where positive is in the direction of blade rotation
%Calculate the same thing for the meshgrid
u_m = u  +  q*z_m-r*y_m;
v_m = v  +  r*x_m-p*z_m - RPM/60*2*pi*z_rh_m; 
w_m = w  +  p*y_m-q*x_m + RPM/60*2*pi*y_rh_m;
n_m = -u_m; %affect of rotor velocity on normal relative water velocity in the x direction
t_m = -w_m.*sin(ang_m) - v_m.*cos(ang_m); %affect of rotor velocity

%******************************************************************
%This is the undisturbed water velocity at the rotor blade and mesh grid locations
%******************************************************************
%calculate the water velocity at the blade and meshgrid locations (10) V_o
%water wave velocity vector
N=length(Z_b(:,1));

if waves.H==0 %if the significant wave height is set to zero then orbital velocities are not calculated
    U_wave=zeros(rb.N_blades,N);
    V_wave=zeros(rb.N_blades,N);
    W_wave=zeros(rb.N_blades,N);
    U_wave_m=zeros(rb.N_mesh,N);
    V_wave_m=zeros(rb.N_mesh,N);
    W_wave_m=zeros(rb.N_mesh,N);
else
   
    X_b2 = [];
    Y_b2 = [];
    Z_b2 = [];
    X_m2 = [];
    Y_m2 = [];
    Z_m2 = [];
    
    for i=1:rb.N_blades %I think that this creates a vector 
      X_b2=[X_b2,X_b(:,i)'];
      Y_b2=[Y_b2,Y_b(:,i)'];
      Z_b2=[Z_b2,Z_b(:,i)'];
    end
    for i=1:N_mesh 
      X_m2=[X_m2,X_m(:,i)'];
      Y_m2=[Y_m2,Y_m(:,i)'];
      Z_m2=[Z_m2,Z_m(:,i)'];
    end

    [teta2,X_vec_2]=meshgrid(waves.teta2,X_b2);
    [k2,Y_vec_2]=meshgrid(waves.k2,Y_b2);
    [w2,Z_vec_2]=meshgrid(waves.w2,Z_b2);
    [H,xx]=meshgrid(waves.H,ones(1,rb.N_blades*N));
    [pri,xx]=meshgrid(waves.pri,ones(1,rb.N_blades*N));
    
    [teta_m2,X_vec_m2]=meshgrid(waves.teta2,X_m2);
    [k_m2,Y_vec_m2]=meshgrid(waves.k2,Y_m2);
    [w_m2,Z_vec_m2]=meshgrid(waves.w2,Z_m2);
    [H_m,xx]=meshgrid(waves.H,ones(1,N_mesh*N));
    [pri_m,xx]=meshgrid(waves.pri,ones(1,N_mesh*N));
    
    UV_wave2=H.*waves.g.*k2./w2.*cosh(k2.*(waves.h-Z_vec_2))./cosh(k2.*waves.h).*sin(k2.*(X_vec_2.*cos(teta2)+Y_vec_2.*sin(teta2))-w2.*t+pri);
    U_wave2=sum(cos(teta2).*UV_wave2,2)'; %calculates the forward velocity of the particle of water
    V_wave2=sum(sin(teta2).*UV_wave2,2)'; %calculates the forward velocity of the particle of water
    W_wave2=sum(H.*waves.g.*k2./w2.*sinh(k2.*(waves.h-Z_vec_2))./cosh(k2.*waves.h).*cos(k2.*(X_vec_2.*cos(teta2)+Y_vec_2.*sin(teta2))-w2.*t+pri),2)'; %calculates the downward velocity of the particle of water

    UV_wave_m2=H_m.*waves.g.*k_m2./w_m2.*cosh(k_m2.*(waves.h-Z_vec_m2))./cosh(k_m2.*waves.h).*sin(k_m2.*(X_vec_m2.*cos(teta_m2)+Y_vec_m2.*sin(teta_m2))-w_m2.*t+pri_m);
    U_wave_m2=sum(cos(teta_m2).*UV_wave_m2,2)'; %calculates the forward velocity of the particle of water
    V_wave_m2=sum(sin(teta_m2).*UV_wave_m2,2)'; %calculates the forward velocity of the particle of water
    W_wave_m2=sum(H_m.*waves.g.*k_m2./w_m2.*sinh(k_m2.*(waves.h-Z_vec_m2))./cosh(k_m2.*waves.h).*cos(k_m2.*(X_vec_m2.*cos(teta_m2)+Y_vec_m2.*sin(teta_m2))-w_m2.*t+pri_m),2)'; %calculates the downward velocity of the particle of water

    
    for i=1:rb.N_blades
        U_wave(i,:)=U_wave2((1+(i-1)*N):(i*N));
        V_wave(i,:)=V_wave2((1+(i-1)*N):(i*N));
        W_wave(i,:)=W_wave2((1+(i-1)*N):(i*N));
    end
    for i=1:N_mesh
        U_wave_m(i,:)=U_wave_m2((1+(i-1)*N):(i*N));
        V_wave_m(i,:)=V_wave_m2((1+(i-1)*N):(i*N));
        W_wave_m(i,:)=W_wave_m2((1+(i-1)*N):(i*N));
    end
end

%% calculate the turbulence over the swept area of the rotor blade here
u_1=sum(turbulence.mag_x(:,2:end-1).*sin(2*pi*(turbulence.fm'*ones(1,turbulence.M-2))*t+turbulence.ang_x(:,2:end-1))); %2 values (1st and last) are removed because these are used in body and VBC
v=sum(turbulence.mag_y(:,2:end-1).*sin(2*pi*(turbulence.fm'*ones(1,turbulence.M-2))*t+turbulence.ang_y(:,2:end-1)));
w=sum(turbulence.mag_z(:,2:end-1).*sin(2*pi*(turbulence.fm'*ones(1,turbulence.M-2))*t+turbulence.ang_z(:,2:end-1)));
u=u_1+turbulence.r_uv*v+turbulence.r_uw*w;

%reshape u,v and w to add with wave velocities
u_adjusted=reshape(u,N_mesh,length(rb.radius))';
%u_adjusted=u_adjusted';
v_adjusted=reshape(v,N_mesh,length(rb.radius))';
%v_adjusted=v_adjusted';
w_adjusted=reshape(w,N_mesh,length(rb.radius))';
%w_adjusted=w_adjusted';


%******TRANSITION FROM MESHGRID TO ROTOR BLADE******
%find index of mesh angles that are just greater than blade angles 
%ang_mesh = 2*pi*(0:1/N_mesh:1); %mesh from 0 to 2pi  
ang_blades2 = ang_global; %angle of the blades
ang_blades1 = atan2(sin(ang_blades2),cos(ang_blades2)); %angle of the blades from -pi to pi
ang_blades = 2*pi*sign(abs(ang_blades1) - ang_blades1) + ang_blades1; %angle of blades from 0 - 2*pi
ang_b_I = ang_blades/(2*pi/N_mesh); %blade angle scaled with respect to mesh grid
index_b_up = ceil(ang_b_I); %Index of mesh grid above blade angle
index_b_down = floor(ang_b_I); %Index of mesh grid below blade angle

%conduct linear interpolation to the blades as they rotate through the mesh
up_sf = ((index_b_up-ang_b_I)'*ones(1,length(r_sec)))'; %normalized angle between the above mesh angle and blade angle
down_sf = ((1+ang_b_I-index_b_up)'*ones(1,length(r_sec)))'; %normalized angle between the below mesh angle and blade angle
long_index_b_down = index_b_down+1;
long_index_b_up = index_b_up+1;
long_vec = [N_mesh 1:N_mesh];
index_b_down_new = long_vec(long_index_b_down);
index_b_up_new = long_vec(long_index_b_up);

u_adjusted_b = u_adjusted(:,index_b_up_new).*down_sf+u_adjusted(:,index_b_down_new).*up_sf;
v_adjusted_b = v_adjusted(:,index_b_up_new).*down_sf+v_adjusted(:,index_b_down_new).*up_sf;
w_adjusted_b = w_adjusted(:,index_b_up_new).*down_sf+w_adjusted(:,index_b_down_new).*up_sf;

ADCP.inflowSelection = 1;
%water total water velocity in the inertial coordinate system
%InflowSelection=1 if using an imported dataset else use uniform values.
if ADCP.inflowSelection == 1 
    time_in = ADCP.tStart+datenum([0,0,0,0,0,t]);
    Uw_interp_b_vec = interp2(ADCP.z,ADCP.t,ADCP.u',Z_b(:),time_in);
    Vw_interp_b_vec = interp2(ADCP.z,ADCP.t,ADCP.v',Z_b(:),time_in);
    Ww_interp_b_vec = interp2(ADCP.z,ADCP.t,ADCP.w',Z_b(:),time_in);
    Uw_interp_m_vec = interp2(ADCP.z,ADCP.t,ADCP.u',Z_m(:),time_in);
    Vw_interp_m_vec = interp2(ADCP.z,ADCP.t,ADCP.v',Z_m(:),time_in);
    Ww_interp_m_vec = interp2(ADCP.z,ADCP.t,ADCP.w',Z_m(:),time_in);
    Uw_interp_b = reshape(Uw_interp_b_vec,size(Z_b)); 
    Vw_interp_b = reshape(Vw_interp_b_vec,size(Z_b)); 
    Ww_interp_b = reshape(Ww_interp_b_vec,size(Z_b));
    Uw_interp_m = reshape(Uw_interp_m_vec,size(Z_m)); 
    Vw_interp_m = reshape(Vw_interp_m_vec,size(Z_m)); 
    Ww_interp_m = reshape(Ww_interp_m_vec,size(Z_m));
    Uw_b = Uw_interp_b+U_wave'+u_adjusted_b;
    Vw_b = Vw_interp_b+V_wave'+v_adjusted_b;
    Ww_b = Ww_interp_b+W_wave'+w_adjusted_b;
    Uw_m = Uw_interp_m+U_wave_m'+u_adjusted;
    Vw_m = Vw_interp_m+V_wave_m'+v_adjusted;
    Ww_m = Ww_interp_m+W_wave_m'+w_adjusted;
else
    Uw_b = U_vec(1)+Z_b*U_vec(4)+U_wave'+u_adjusted_b; %North Water Velocity -replace first 2 terms with U_adcp
    Vw_b = U_vec(2)+Z_b*U_vec(5)+V_wave'+v_adjusted_b; %East Water Velocity
    Ww_b = U_vec(3)+Z_b*U_vec(6)+W_wave'+w_adjusted_b; %Down Water Velocity
    Uw_m = U_vec(1)+Z_m*U_vec(4)+U_wave_m'+u_adjusted;
    Vw_m = U_vec(2)+Z_m*U_vec(5)+V_wave_m'+v_adjusted;
    Ww_m = U_vec(3)+Z_m*U_vec(6)+W_wave_m'+w_adjusted;
end

u_wv_b = (L_BI(1,1)*Uw_b' + L_BI(1,2)*Vw_b' + L_BI(1,3)*Ww_b')'; %forward water velocity (undisturbed)
v_wv_b = (L_BI(2,1)*Uw_b' + L_BI(2,2)*Vw_b' + L_BI(2,3)*Ww_b')';
w_wv_b = (L_BI(3,1)*Uw_b' + L_BI(3,2)*Vw_b' + L_BI(3,3)*Ww_b')';
u_wv_m = (L_BI(1,1)*Uw_m' + L_BI(1,2)*Vw_m' + L_BI(1,3)*Ww_m')';
v_wv_m = (L_BI(2,1)*Uw_m' + L_BI(2,2)*Vw_m' + L_BI(2,3)*Ww_m')';
w_wv_m = (L_BI(3,1)*Uw_m' + L_BI(3,2)*Vw_m' + L_BI(3,3)*Ww_m')';

%sp_wv_b = sqrt(u_wv_b.^2+u_wv_b.^2+u_wv_b.^2);
n_w_b = u_wv_b; %normal undisturbed water velocity (is defined positive in the x direction)
t_w_b = w_wv_b.*sin(ang_b)+v_wv_b.*cos(ang_b); %tangental undisturbed water velocity (toward leadind edge of blade)
sp_wv_m = sqrt(u_wv_m.^2+u_wv_m.^2+u_wv_m.^2);
n_w_m = u_wv_m; %normal undisturbed water velocity (is defined positive in the x direction)
t_w_m = w_wv_m.*sin(ang_m)+v_wv_m.*cos(ang_m); %tangental undisturbed water velocity (toward leadind edge of blade)


%**************************************************************************
%this it the total water velocity and forces
%**************************************************************************
%normal and tangental actual water velocities

n_rel_b = n_w_b + n_b + W_n_b_old; %normal water + (-rotor vel) + (-wake reduction?)
t_rel_b = t_w_b + t_b + W_t_b_old;
n_rel_m = n_w_m + n_m + W_n_m_old;
t_rel_m = t_w_m + t_m + W_t_m_old;

%calculate the angles of attack
ang_w = atan2(-n_rel_b,-t_rel_b);
ang_w_m = atan2(-n_rel_m,-t_rel_m);

%b_b is radians for 4pi/3+ang 2pi/3+ang ang
att = atan2(sin(ang_w - b_b),cos(ang_w - b_b));
att_m = atan2(sin(ang_w_m - b_m),cos(ang_w_m - b_m));
%calculate the coeficients of lift and drag for the rotor blade sections
%from coefficient matrix / at discrete points 

angles_mat = pi/180*(rb.angles*ones(1,25));
element_mat = ones(length(rb.angles),1)*(1:25);
element_att = (1:25)'*ones(1,3);
element_att_m = (1:25)'*ones(1,N_mesh);
cl = interp2(element_mat,angles_mat,rb.Cl_matrix,element_att,att);
cd = interp2(element_mat,angles_mat,rb.Cd_matrix,element_att,att);
cl_m = interp2(element_mat,angles_mat,rb.Cl_matrix,element_att_m,att_m);


%calculate the lift and drag off of each section (16)
L = 0.5*mdl.rho*((rb.rad_thick.*rb.chord)*ones(1,3)).*cl.*(n_rel_b.^2+t_rel_b.^2);
D = 0.5*mdl.rho*((rb.rad_thick.*rb.chord)*ones(1,3)).*cd.*(n_rel_b.^2+t_rel_b.^2);

L_m = 0.5*mdl.rho*((rb.rad_thick.*rb.chord)*ones(1,N_mesh)).*cl_m.*(n_rel_m.^2+t_rel_m.^2);

%convert lift and drag to axial (normal) and tangental loads (18)
A = D.*sin(ang_w)+L.*cos(ang_w);
T = -D.*cos(ang_w)+L.*sin(ang_w);

%convert the forces to the body fixed frame (18)
fx = -A;
fy = T.*cos(ang_b);
fz = T.*sin(ang_b);

%analyze force on a section of a blade due to turbulence
fx_oneblade=sum(fx(:,1));

%these are outputs for the rotor analysis
MA_blade1 = sum(rb.radius.*A(:,1));
MA_blade2 = sum(rb.radius.*A(:,2));
MA_blade3 = sum(rb.radius.*A(:,3));
%MT_blade1 = sum(rb.radius.*T(:,1));
    
%calculate the total forces in the body fixed frame and moments on the
%rotor about the shaft (19)
fx_t = sum(sum(fx'));
fy_t = sum(sum(fy'));
fz_t = sum(sum(fz'));

%calculate the total moments on the rotor about the shaft (21)
Mx_s = sum(sum(fz.*y_rh-fy.*z_rh));%-sum(A1.*r_sec+A2.*r_sec+A3.*r_sec);
My_s = sum(sum(fx.*z_rh-fz.*0));
Mz_s = sum(sum(fy.*0-fx.*y_rh));

%calculate the total moments on the turbine about its CG (20)
Mx = sum(sum(fz.*y_b-fy.*z_b));
My = sum(sum(fx.*z_b-fz.*x_b));
Mz = sum(sum(fy.*x_b-fx.*y_b));

%**************************************************************************
%this is the updated induced velocity vector (ONLY FOR MESHGRID)
%**************************************************************************
%Prandtl's tip loss factor
%calculate Prandtl's tip loss correlation factor (13)
R = 10;
F = 2/pi*acos(exp(-(rb.N_blades*(R-r_m)./(2*r_m.*sin(abs(ang_w_m))+eps))));

%axial induction factor
a = W_0_n_m_old./sqrt(Uw_m.^2 + Vw_m.^2 + Ww_m.^2+eps);%this should be positive

%Glauert correction factor
a_c = 0.2;
fg_a = sign((a_c-a)+abs(a_c-a)); %this is for a_c > a
fg = (1-fg_a).*((a_c./(a+eps)).*(2-a_c./(a+eps))) + fg_a;%%%%%%%%%%%%%%%%%%%

%quasi steady curent generated by the turbine
L_ul = 0.5*mdl.rho*(rb.chord*ones(1,N_mesh)).*cl_m.*(n_rel_m.^2+t_rel_m.^2);
W_qs_n_m = rb.N_blades*L_ul.*cos(ang_w_m)./(4*pi*mdl.rho*r_m.*F.*sqrt((u_wv_m+fg.*W_0_n_m_old).^2+v_wv_m.^2+w_wv_m.^2)+eps);
W_qs_t_m = -rb.N_blades*L_ul.*sin(ang_w_m)./(4*pi*mdl.rho*r_m.*F.*sqrt((u_wv_m+fg.*W_0_n_m_old).^2+v_wv_m.^2+w_wv_m.^2)+eps);

%H factor for unsteady induced velocities
k = 0.6;
%a should not be allowed to exceed 0.5
overp5 = 0.5*sign((a-0.5)+abs(a-0.5));
underp5 = (1-2*overp5).*a;
a_limited = (overp5 + underp5).*(sign(overp5 + underp5)+1)/2;
tau1 = (1.1./(1-1.3*a_limited)).*(R./(sp_wv_m+eps));
H_n = W_qs_n_m + k*tau1.*(W_qs_n_m-W_qs_n_m_old)/(delta_t+eps);
H_t = W_qs_t_m + k*tau1.*(W_qs_t_m-W_qs_t_m_old)/(delta_t+eps);

%intermediate W value
W_int_n_m = H_n + (W_int_n_m_old-H_n).*exp(-delta_t./tau1);
W_int_t_m = H_t + (W_int_t_m_old-H_t).*exp(-delta_t./tau1);

%new value of W without the yaw factored in
tau2 = (0.39-0.26*(r_m/R).^2).*tau1;
W_0_n_m = W_int_n_m + (W_0_n_m_old - W_int_n_m).*exp(-delta_t./tau2);
W_0_t_m = W_int_t_m + (W_0_t_m_old - W_int_t_m).*exp(-delta_t./tau2);

%calculate the new W with the yaw factored in
W_0_x_m = W_0_n_m;
W_0_y_m = W_0_t_m.*cos(ang_m);
W_0_z_m = W_0_t_m.*sin(ang_m);

%******TRANSITION FROM MESHGRID TO ROTOR BLADE******
%find index of mesh angles that are just greater than blade angles 
%ang_mesh = 2*pi*(0:1/N_mesh:1); %mesh from 0 to 2pi  
% ang_blades2 = ang_global; %angle of the blades
% ang_blades1 = atan2(sin(ang_blades2),cos(ang_blades2)); %angle of the blades from -pi to pi
% ang_blades = 2*pi*sign(abs(ang_blades1) - ang_blades1) + ang_blades1; %angle of blades from 0 - 2*pi
% ang_b_I = ang_blades/(2*pi/N_mesh); %blade angle scalesd with respect to mesh grid
% index_b_up = ceil(ang_b_I); %Index of mesh grid above blade angle
% index_b_down = floor(ang_b_I); %Index of mesh grid below blade angle
% 
% %conduct linear interpolation
% up_sf = ((index_b_up-ang_b_I)'*ones(1,length(r_sec)))'; %normalized angle between the above mesh angle and blade angle
% down_sf = ((1+ang_b_I-index_b_up)'*ones(1,length(r_sec)))'; %normalized angle between the below mesh angle and blade angle
% %W_0_x_b = W_0_x_m(:,index_b_up).*down_sf+W_0_x_m(:,index_b_down).*up_sf;
% %W_0_y_b = W_0_y_m(:,index_b_up).*down_sf+W_0_y_m(:,index_b_down).*up_sf;
% %W_0_z_b = W_0_z_m(:,index_b_up).*down_sf+W_0_z_m(:,index_b_down).*up_sf;
% long_index_b_down = index_b_down+1;
% long_index_b_up = index_b_up+1;
% long_vec = [N_mesh 1:N_mesh];
% index_b_down_new = long_vec(long_index_b_down);
% index_b_up_new = long_vec(long_index_b_up);

W_0_n_b = W_0_n_m(:,index_b_up_new).*down_sf+W_0_n_m(:,index_b_down_new).*up_sf;
W_0_t_b = W_0_t_m(:,index_b_up_new).*down_sf+W_0_t_m(:,index_b_down_new).*up_sf;


%******THIS IS FOR BOTH THE ROTOR BLADE AND MESH****************
%This is the wake angle skew angle defined as the angle between the water
%velocity in the wake at r/R = 0.7 averaged over the three blades
[val, index] = min(abs(r_sec/R-0.7));
%ang_sk = atan2(sqrt((sum(W_0_y(index,:)+v_wv_rb(index,:)))^2+(sum(W_0_z(index,:)+w_wv_rb(index,:)))^2),-sum((W_0_x(index,:)+u_wv_rb(index,:))));
ang_sk = atan2(sqrt((sum(W_0_y_m(index,:)+v_wv_m(index,:)))^2+(sum(W_0_z_m(index,:)+w_wv_m(index,:)))^2),-sum((W_0_x_m(index,:)+u_wv_m(index,:))));
%this is the angle of the deepest portion of the rotor
%y and z velocity of the water with respect to the shaft at the shaft

v_rel_mat = t_rel_m.*cos(ang_m);
v_rel = sum(v_rel_mat(index,:))/N_mesh;
w_rel_mat = t_rel_m.*sin(ang_m);
w_rel = sum(w_rel_mat(index,:))/N_mesh;
ang_deep = atan2(v_rel,-w_rel);

W_n_m = W_0_n_m.*(1+r_m/R.*tan(ang_sk/2).*cos(ang_m-ang_deep));
W_t_m = W_0_t_m.*(1+r_m/R.*tan(ang_sk/2).*cos(ang_m-ang_deep));

W_n_b = W_0_n_b.*(1+r_b/R.*tan(ang_sk/2).*cos(ang_b-ang_deep));
W_t_b = W_0_t_b.*(1+r_b/R.*tan(ang_sk/2).*cos(ang_b-ang_deep));
%f_blades = sqrt((sum(fx)).^2+(sum(fy)).^2+(sum(fz)).^2);

rt_m_g1 = root_twist_m_global(1); 
rt_m_g5 = root_twist_m_global(5);