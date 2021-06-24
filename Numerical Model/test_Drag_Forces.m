
%% Below here are constants from user inputs
N_mesh=20;
L_cable=1;
D_cable=1;
m_cable=1;
disp_m_cable = 1;
EA_cable = 1;
Cd_cable = 1;
U = 1.5;
dUdZ = 0;
N_cable = 5;
Hs = 0;
TI = 0.001;
Z_disp = 0;
Y_disp = 0;

%% Below here is from the mask

[mdl rb] = turbine_constants_20mD_turbulence_VarBuoy_v2(N_mesh);
Z0_ZCAP = 325;
%initial try for v11
%this is for 1.6 m/s
state0 = [0, 0, 0, 0, 0, 0, 555.8476,1.6924+Y_disp, 50.3603+Z_disp, 0.0073, 0.0189, 3.1485];

Current=[U; 0; 0; dUdZ; 0; 0];%Current speed vector

first_node = [0 0 0 0 0 Z0_ZCAP];
last_node = attachment_point([state0(1:4),0,state0(5:10),0,state0(11:12)],mdl);
if N_cable == 1;
    state_cable0 = []; %intermediate nodes
elseif N_cable == 2;
    state_cable0 = 0.5*last_node+0.5*first_node;
elseif N_cable == 3;
    state_cable0 = [0.3333*last_node+0.6666*first_node 0.6666*last_node+0.3333*first_node];
elseif N_cable == 4;
    state_cable0 = [0.25*last_node+0.75*first_node 0.5*last_node+0.5*first_node 0.75*last_node+0.25*first_node];
elseif N_cable == 5;
    state_cable0 = [0.2*last_node+0.8*first_node 0.4*last_node+0.6*first_node 0.6*last_node+0.4*first_node 0.8*last_node+0.2*first_node];
elseif N_cable == 8;
    state_cable0 = [0.125*last_node+0.8750*first_node 0.25*last_node+0.75*first_node 0.375*last_node+0.625*first_node 0.5*last_node+0.5*first_node 0.625*last_node+0.375*first_node 0.75*last_node+0.25*first_node 0.878*last_node+0.125*first_node];
else
    disp('cable model curently only valid for 1-5 elements')
end
    
Mod_vars.l = L_cable;
Mod_vars.D = D_cable;
Mod_vars.m = m_cable;
Mod_vars.disp_m = disp_m_cable;
Mod_vars.N_cable = N_cable;
Mod_vars.EA = EA_cable;
Mod_vars.Cd = Cd_cable;
Mod_vars.drag_t = 1/4*mdl.rho*2*pi*Mod_vars.D/2*Mod_vars.l*Mod_vars.Cd;
Mod_vars.drag_n = 1/4*mdl.rho*Mod_vars.D*Mod_vars.l*Mod_vars.Cd;
Mod_vars.K = Mod_vars.EA/Mod_vars.l;
Mod_vars.Wt = (Mod_vars.m*mdl.g-Mod_vars.disp_m*mdl.g)*Mod_vars.l;
Mod_vars.M = Mod_vars.m*Mod_vars.l;
Mod_vars.Cid = 1000;
Mod_vars.rho = mdl.rho;
Mod_vars.Z0_ZCAP = Z0_ZCAP;
Mod_vars.Current=Current;

%% This calls the wave spectrum that is used
n=20; %number of wave elements
m=5;  %number of wave directions
%Hs=0;%1.859; %significant wave height
angle0=0; %mean direction of wave propagation
s=25; %not sure
waves = waves_spec_calc2(n,m,Hs,angle0,s);

%% This is what calculates all of the turbulence constants
%TI = 0.10; %Turbulence intensity
a=1;b=1; %sigma_y=a*sigma_x,sigma_z=b*sigma_x
fmin = 0.01;%minimum frequency to resolve
fmax=1; %maximum frequency to resolve
C=27; %coherence decay constant, Veers (1986)
% N=round(2*fmax*T*1.2); %data points
P=1500;%number of frequency discretization
turbulence = turbulence_constants_v3(TI,fmin,fmax,rb.radius,N_mesh,U,a,b,C,P); %this creates the required wave data

%% This section processes ADCP data and outputs current velocities
filename = 'ADCP-Data-Irma.mat'; %Input file with ADCP data
tStart = datenum('15-Sep-0017 04:32:00');%start time for generating a subset of data
tStop = datenum('03-Nov-0017 18:14:00');%stop time for generating a subset of data
ADCP = ADCP_Init(filename,tStart,tStop)

%% This provides the Simulink inputs and calls the hydrodynamics functions
% u = in(1);
% v = in(2);
% w = in(3);
% p = in(4);
% p_r = in(5);
% q = in(6);
% r = in(7);
% X = in(8);
% Y = in(9);
% Z = in(10);
% phi = in(11);
% phi_r = in(12);
% theta = in(13);
% psi = in(14);
% Uw = in(15); %North Water Velocity
% Vw = in(16); %East Water Velocity
% Ww = in(17); %Down Water Velocity
% dUw = in(18); %North Water Velocity Gradient
% dVw = in(19); %East Water Velocity Gradient
% dWw = in(20); %Down Water Velocity Gradient
% t=in(21);
% t_delayed = in(22);
% B3_ang_rad = in(23);
% B2_ang_rad = in(24);
% B1_ang_rad = in(25);
in = [0 0 0 0 1 0 0 0 0 60 0 0 0 0 1 0 0 0 0 0 0.1 0 0 0 0] 
[out] = drag_forces_20m_IBP_v2(in,mdl,rb,waves,turbulence,ADCP)
