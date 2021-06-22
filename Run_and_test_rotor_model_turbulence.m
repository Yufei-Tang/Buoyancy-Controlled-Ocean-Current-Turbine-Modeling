clear all
close all
clc

%% calculate mdl and rb
% these are the constant system model (mdl) and rotor blade model (rb)
% properties
N_mesh = 32;
%this loads all of the constants that define the turbine
[mdl rb] = turbine_constants_20mD_turbulence(N_mesh);
%associated sub-programs and other files include
%matrix = xlsread('TURB_0000_Hydrostatic_Model_110308.xls','Turbine Components','O3:AB40');
%rb = rotor_constants_20mD_v1(N_mesh);


%% These are the parameters that define the wave field (waves)
n=10; %this is the number of sin waves used for each wave angle (around 25 is good)
m=3; %this is the number of propication directions for the waves
Hs=0; %SIGNIFICANT WAVE HEIGHT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
angle=0; %mean angle a wave propication (0 means the waves are headed north)
s=25; %spreading factor for the wave fild (should be around 25)
waves = waves_spec_calc2(n,m,Hs,angle,s); %this creates the required wave data


%% These are the states of the turbine
%This is the yaw angle of the rotor (we should probably keep this 0 for most of our analyses)
yaw = 0;
%This is the flow speed in meters/second !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
U = 1.6;
%THIS IS THE DEPTH OF THE CENTER OF THE ROTOR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Z = 20;
turb_states = [0 0 0 0 0 0 0 0 Z 0 0 yaw];
water_vel = [-U 0 0 0 0 0];
%root_twist = [0 0 0]*pi/180;


%% This is what calculates all of the turbulence constants
%put the user constants here
TI = 0.20; %Turbulence intencity
a=1;b=1; %sigma_y=a*sigma_x,sigma_z=b*sigma_x
fmin = 0.01;%minimum frequency to resolve
fmax=1; %maximum frequency to resolve
C=5; %coherence decay constant, Veers (1986)
% N=round(2*fmax*T*1.2); %data points
P=1000;%number of frequency discretization
turbulence = turbulence_constants_v2(TI,fmin,fmax,rb.radius,N_mesh,U,a,b,C,P); %this creates the required turbulence data


%% This runs the numerical rotor model
TSR_max = 10.14;
RPM = TSR_max*U*60/2/pi/10; %this is the RPM for max TSR
total_time = 60*5+30; % total time for the run
time_step = 0.05; % time step
ct = 0;
time_old = 0;
for time = time_step:time_step:total_time
    ct = ct + 1;
    ang = 2*pi/60*RPM*time;  % rotor rotation angle about axis in radians (this will make the rotor angles vary appropriately for a fixed RPM)  
%     [fx_t, fy_t, fz_t, Mx, My, Mz, Mx_s, My_s, Mz_s, MA_blade1, MA_blade2, MA_blade3] = rotor_forces_20mD_turbulence([turb_states time time_old],water_vel,[RPM ang],mdl,rb,waves,turbulence,[0 0 0]);
    [fx_oneblade fx_t, fy_t, fz_t, Mx, My, Mz, Mx_s, My_s, Mz_s, MA_blade1, MA_blade2, MA_blade3]= rotor_forces_20mD_turbulence([turb_states time time_old],water_vel,[RPM ang],mdl,rb,waves,turbulence,[0 0 0]);
    time_old = time;

    %This just saves the variables of interest as vectors
    ang_vec(ct) = ang;
    My_vec(ct) = My; 
    shaft_torque(ct) = Mx_s;
    shaft_power(ct) = RPM*2*pi/60*Mx_s;
    rotor_drag(ct) = fx_t;
    oneblade_drag(ct)=fx_oneblade;
    time_vector(ct) = time;
end




% figure(1)
% plot(time_vector-30,shaft_torque/1000)
% axis([0 30 200 250])
% xlabel('time (s)')
% ylabel('shaft torque (kN*m)')
% 
% figure(2)
% plot(time_vector-30,shaft_power/1000)
% axis([0 30 300 400])
% xlabel('time (s)')
% ylabel('shaft power (kW)')

figure(3)
plot(time_vector-30,shaft_torque/1000,'r',time_vector-30,shaft_power/1000,':b')
legend('shaft torque','shaft power')
axis([0 total_time-30 175 400])
xlabel('time (s)')
ylabel('shaft power (kW)/shaft torque (kN-m)')
title('TI=0.10')

avg_power=mean(shaft_power(600:length(time_vector)))
avg_torque=mean(shaft_torque(600:length(time_vector)))

figure(4)
plot(time_vector-30, -oneblade_drag/1000)
axis([0 total_time-30 120 140])
xlabel('time (s)')
ylabel('force on one blade (kN)')

% figure(5)
% plot(time_vector,-rotor_drag/1000)
% xlabel('time (s)')
% ylabel('Drag Force (-fx) (kN)')


figure(3)
plot(tout-30,Torque/1000,'r')
legend('shaft torque')
axis([0 tout-30 175 140])
xlabel('time (s)')
ylabel('shaft power (kW)/shaft torque (kN-m)')
title('TI=0.10')

figure(4)
plot(tout, fx_oneblade/1000)
axis([0 tout-30 120 140])
xlabel('time (s)')
ylabel('force on one blade (kN)')
