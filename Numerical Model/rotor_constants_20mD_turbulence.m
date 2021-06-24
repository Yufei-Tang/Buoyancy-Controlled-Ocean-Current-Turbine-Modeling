function rb = rotor_constants_20mD_turbulence(N_mesh)

%% this is the mass of the rotor scaled from the Verdant Rotor 
mass = (415/2.2)*(20/3)^3; %mass of the rotor in kg
rho_steel = 7850; %density of steel

%% this is the moment of inertia for the rotor hub  
N_blades = 3; %this is the number of blades (it is not used consistantly though the code yet)
D_ring = 22.5*0.0254;%diamiter of the inner ring in m
L_ring = 12.5*0.0254;%length or depth of the inner ring
t_ring = 3/8*0.0254;
m_ring = L_ring*t_ring*D_ring*pi*rho_steel;
Ixx_ring = m_ring*D_ring^2/4;

N_plates = 3;
ID_plates = 9*0.0254;
OD_plates = 22*0.0254;
t_plates = 3/8*0.0254;
fraction_plates = 0.5; %frontal area fraction
m_plates = N_plates*(OD_plates^2-ID_plates^2)/4*pi*t_plates*fraction_plates*rho_steel;
Ixx_plates = 0.5*m_plates*(OD_plates^2+ID_plates^2)/4;

m_hub = (m_ring+m_plates)*(20/3)^3; %20 m
Ixx_hub = (Ixx_ring+Ixx_plates)*(20/3)^5; %20 m

m_blades = mass-m_hub; %20 m

load('rotor_data_for_20m_fx83_varSp_varPi_P2F_GS.mat')
%this calculates the radual thickness of each section
rad_thick(1) = (rb.radius(2)-rb.radius(1))/2+(rb.radius(1)-1.5);%20 m
rad_thick(2:24) = (rb.radius(3:25)-rb.radius(1:23))/2; %20 m
rad_thick(25) = (rb.radius(25)-rb.radius(25))/2+(10-rb.radius(25)); %20 m
rb.rad_thick = rad_thick'; %20 m
m_sections = (rb.rad_thick.*rb.chord.*rb.thick)/sum(rb.rad_thick.*rb.chord.*rb.thick)*m_blades; %20 m
Ixx_blade = sum(m_sections.*rb.radius.^2); %20 m

%these vectors are output for the simulation
rb.N_blades=N_blades;
rb.N_mesh = N_mesh;
rb.Ixx = Ixx_hub+Ixx_blade; %20 m

%global W_n_b_old W_t_b_old W_n_m_old W_t_m_old W_0_n_m_old W_0_t_m_old W_int_n_m_old W_int_t_m_old W_qs_n_m_old W_qs_t_m_old;
%global variables used in the rotor model
global W_qs_n_m W_qs_t_m W_int_n_m W_int_t_m W_0_n_m W_0_t_m W_n_m W_t_m W_n_b W_t_b time_delayed_old root_twist_m_global ang_global
W_qs_n_m = zeros(length(rb.radius),N_mesh);  %quasi-steady wake fraction in the normal direction
W_qs_t_m = zeros(length(rb.radius),N_mesh);  %quasi-steady wake fraction in the tangental direction
W_int_n_m = zeros(length(rb.radius),N_mesh); %intermediate wake fraction in the normal direction
W_int_t_m = zeros(length(rb.radius),N_mesh); %intermediate wake fraction in the tangental direction
W_0_n_m = zeros(length(rb.radius),N_mesh);   %wake fraction not compensated for non-axial flows in the normal direction
W_0_t_m = zeros(length(rb.radius),N_mesh);   %wake fraction not compensated for non-axial flows in the tangental direction
W_n_m = zeros(length(rb.radius),N_mesh);     %wake fraction compensated for non axial flows in the normal direction
W_t_m = zeros(length(rb.radius),N_mesh);     %wake fraction compensated for non axial flows in the tangental direction

W_n_b = zeros(length(rb.radius),N_blades);   %wake fraction compensated for non axial flows in the normal direction
W_t_b = zeros(length(rb.radius),N_blades);   %wake fraction compensated for non axial flows in the tangental direction
time_delayed_old = 0;

root_twist_m_global = zeros(1,N_mesh);       %root twist angle for the momentum mesh blades before a rotor blade passes them
ang_global = 2*pi*(1/rb.N_blades:1/rb.N_blades:1); %this is the angle of the tree blades 

global W_qs_n_m_old W_qs_t_m_old W_int_n_m_old W_int_t_m_old W_0_n_m_old W_0_t_m_old W_n_m_old W_t_m_old W_n_b_old W_t_b_old root_twist_m_global_old ang_global_old
W_qs_n_m_old = W_qs_n_m;
W_qs_t_m_old = W_qs_t_m;
W_int_n_m_old = W_int_n_m;
W_int_t_m_old = W_int_t_m;
W_0_n_m_old = W_0_n_m;
W_0_t_m_old = W_0_t_m;
W_n_m_old = W_n_m;
W_t_m_old = W_t_m;
W_n_b_old = W_n_b;
W_t_b_old = W_t_b;
root_twist_m_global_old = zeros(1,N_mesh);       %root twist angle for the momentum mesh blades before a rotor blade passes them
ang_global_old = 2*pi*(1/rb.N_blades:1/rb.N_blades:1); %this is the angle of the three blades
