function [mdl, rb] = turbine_constants_20mD_turbulence_VarBuoy_v2(N_mesh)
mdl.g = 9.81;
mdl.rho = 1025.18;

matrix = xlsread('TURB_0000_Hydrostatic_Model_110308.xls','Turbine Components','O3:AB40');

%this calls the rotor constants code
rb = rotor_constants_20mD_turbulence(N_mesh);

%rotor moments of inertia about shaft force location
Ixx_rotor = rb.Ixx;% + 0; %20 m
Iyy_rotor = 0.5*rb.Ixx;% + mdl.mass_rotor*(mdl.X_rotor_wrt_shaft)^2;%20 m
Izz_rotor = 0.5*rb.Ixx;% + mdl.mass_rotor*(mdl.X_rotor_wrt_shaft)^2;%20 m

%define the origin of the body fixed coordinate system
x_origin = (matrix(25,3)+19.7)*0.0254;
y_origin = matrix(25,4)*0.0254;
z_origin = matrix(25,5)*0.0254;

%change the buoyancy and CB/CG of the Variable Buoyancy Tank
matrix(26,2) = matrix(26,2)-32; %decrease the displaced mass of the BCM by 12.5
matrix(26,3) = matrix(26,3)+7.5; %***************************************** move BCMs 7.5 inches forward on 3 m model
matrix(26,9) = matrix(26,9)+0.9+7.5;%************************************** move BCMs 7.5 inches forward on 3 m model
matrix(26,5) = matrix(26,5)-4.5*(3/20)/0.0254;%**************************** move BCMs 4.5 meters up on 20 m model
matrix(26,11) = matrix(26,11)-4.5*(3/20)/0.0254;%************************** move BCMs 4.5 meters up on 20 m model

%this finds the center of gravity and center of buoyancy of all components
mdl.mass = (matrix(:,1)/2.2)*(20/3)^3; %20 m
mdl.disp_mass = (matrix(:,2)/2.2)*(20/3)^3; %20 m
mdl.Xcg = (matrix(:,3)*0.0254-x_origin)*(20/3);%WRT the Fb coordinate system
mdl.Ycg = (matrix(:,4)*0.0254-y_origin)*(20/3); %20 m
mdl.Zcg = (matrix(:,5)*0.0254-z_origin)*(20/3); %20 m
mdl.Xcb = (matrix(:,9)*0.0254-x_origin)*(20/3); %20 m
mdl.Ycb = (matrix(:,10)*0.0254-y_origin)*(20/3);%20 m
mdl.Zcb = (matrix(:,11)*0.0254-z_origin)*(20/3);%20 m
disp('Variable Buoyancy Chamber: Dry Mass (kg), Disp Mass Fixed Buoyancy (kg), X_cg, Z_cg')
disp(mdl.mass(26))
disp(mdl.disp_mass(26))
disp(mdl.Xcg(26))
disp(mdl.Zcg(26))

%********************ENTIRE SYSTEM*****************************************
%total mass and displaced mass
mass_T = sum(mdl.mass); %20 m
disp_mass_T = sum(mdl.disp_mass); %20 m
% calculate the center of mass and buoyancy for the entire system WRt Fb
mdl.CMx = sum(mdl.Xcg.*mdl.mass)/mass_T; %20 m
mdl.CMy = sum(mdl.Ycg.*mdl.mass)/mass_T; %20 m
mdl.CMz = sum(mdl.Zcg.*mdl.mass)/mass_T; %20 m
mdl.CBx = sum(mdl.Xcb.*mdl.disp_mass)/disp_mass_T; %20 m
mdl.CBy = sum(mdl.Ycb.*mdl.disp_mass)/disp_mass_T; %20 m
mdl.CBz = sum(mdl.Zcb.*mdl.disp_mass)/disp_mass_T; %20 m
% calculate the moments of inertia as if everything is a point mass for the
% entire system pluss the rotor moments of inertia
mdl.Ixx = (sum(mdl.mass.*(mdl.Ycg.^2+mdl.Zcg.^2))+20*(20/3)^5)+Ixx_rotor; %20 m
mdl.Iyy = (sum(mdl.mass.*(mdl.Xcg.^2+mdl.Zcg.^2))+20*(20/3)^5)+Iyy_rotor; %20 m
mdl.Izz = (sum(mdl.mass.*(mdl.Ycg.^2+mdl.Xcg.^2))+20*(20/3)^5)+Izz_rotor; %20 m
%mdl.Ixy = sum(mdl.mass.*(mdl.Ycg.*mdl.Xcg));
mdl.Ixz = sum(mdl.mass.*(mdl.Zcg.*mdl.Xcg)); %20 m
%mdl.Iyz = sum(mdl.mass.*(mdl.Ycg.*mdl.Zcg));
mdl.mass_T = mass_T; %20 m
mdl.disp_mass_T = disp_mass_T; %20 m

%********************MAIN BODY ONLY****************************************
indx_b = [1:23,26:30,33:38];
%total mass and displaced mass
mass_b = sum(mdl.mass(indx_b));  %20 m
disp_mass_b = sum(mdl.disp_mass(indx_b));  %20 m
% calculate the center of mass
mdl.CMx_b = sum((mdl.Xcg(indx_b)).*(mdl.mass(indx_b)))/mass_b; %20 m
mdl.CMy_b = sum((mdl.Ycg(indx_b)).*(mdl.mass(indx_b)))/mass_b; %20 m
mdl.CMz_b = sum((mdl.Zcg(indx_b)).*(mdl.mass(indx_b)))/mass_b; %20 m
mdl.CBx_b = sum((mdl.Xcb(indx_b)).*(mdl.disp_mass(indx_b)))/disp_mass_b; %20 m
mdl.CBy_b = sum((mdl.Ycb(indx_b)).*(mdl.disp_mass(indx_b)))/disp_mass_b; %20 m
mdl.CBz_b = sum((mdl.Zcb(indx_b)).*(mdl.disp_mass(indx_b)))/disp_mass_b; %20 m
% calculate the moments of inertia as if everything is a point mass
mdl.Ixx_b = sum((mdl.mass(indx_b)).*((mdl.Ycg(indx_b)).^2+(mdl.Zcg(indx_b)).^2))+15*(20/3)^5; %20 m
mdl.Iyy_b = sum((mdl.mass(indx_b)).*((mdl.Xcg(indx_b)).^2+(mdl.Zcg(indx_b)).^2))+15*(20/3)^5; %20 m
mdl.Izz_b = sum((mdl.mass(indx_b)).*((mdl.Ycg(indx_b)).^2+(mdl.Xcg(indx_b)).^2))+15*(20/3)^5; %20 m
%mdl.Ixy_b = sum(mdl.mass(indx_b).*(mdl.Ycg(indx_b).*mdl.Xcg(indx_b)));
mdl.Ixz_b = sum((mdl.mass(indx_b)).*((mdl.Zcg(indx_b)).*(mdl.Xcg(indx_b)))); %20 m
%mdl.Iyz_b = sum(mdl.mass(indx_b).*(mdl.Ycg(indx_b).*mdl.Zcg(indx_b)));
mdl.mass_b = mass_b;
mdl.disp_mass_b = disp_mass_b;

%********************ROTOR SECTION ONLY************************************
indx_r = [24,25,31,32];
%total mass and displaced mass
mass_r = sum(mdl.mass(indx_r)); %20 m
disp_mass_r = sum(mdl.disp_mass(indx_r)); %20 m
% calculate the center of mass
mdl.CMx_r = sum(mdl.Xcg(indx_r).*mdl.mass(indx_r))/mass_r; %20 m
mdl.CMy_r = sum(mdl.Ycg(indx_r).*mdl.mass(indx_r))/mass_r; %20 m
mdl.CMz_r = sum(mdl.Zcg(indx_r).*mdl.mass(indx_r))/mass_r; %20 m
mdl.CBx_r = sum(mdl.Xcb(indx_r).*mdl.disp_mass(indx_r))/disp_mass_r; %20 m
mdl.CBy_r = sum(mdl.Ycb(indx_r).*mdl.disp_mass(indx_r))/disp_mass_r; %20 m
mdl.CBz_r = sum(mdl.Zcb(indx_r).*mdl.disp_mass(indx_r))/disp_mass_r; %20 m
% calculate the moments of inertia as if everything is a point mass plus
% the rotor's moment of inertia
mdl.Ixx_r = sum(mdl.mass(indx_r).*(mdl.Ycg(indx_r).^2+mdl.Zcg(indx_r).^2))+5*(20/3)^5+Ixx_rotor; %20 m
mdl.Iyy_r = sum(mdl.mass(indx_r).*(mdl.Xcg(indx_r).^2+mdl.Zcg(indx_r).^2))+5*(20/3)^5+Iyy_rotor; %20 m
mdl.Izz_r = sum(mdl.mass(indx_r).*(mdl.Ycg(indx_r).^2+mdl.Xcg(indx_r).^2))+5*(20/3)^5+Izz_rotor; %20 m
%mdl.Ixy_r = sum(mdl.mass(indx_r).*(mdl.Ycg(indx_r).*mdl.Xcg(indx_r)));
%mdl.Ixz_r = sum(mdl.mass(indx_r).*(mdl.Zcg(indx_r).*mdl.Xcg(indx_r)));
%mdl.Iyz_r = sum(mdl.mass(indx_r).*(mdl.Ycg(indx_r).*mdl.Zcg(indx_r)));
mdl.mass_r = mass_r;  %20 m
mdl.disp_mass_r = disp_mass_r;  %20 m

%***********************INVERSE MATRIX*************************************
mdl.mat_inv = inv(2*...
    [mdl.mass_T, 0, 0, 0, mdl.mass_b*mdl.CMz_b, 0;...
     0, mdl.mass_T, 0, -mdl.mass_b*mdl.CMz_b, 0, mdl.mass_T*mdl.CMx;...
     0, 0, mdl.mass_T, 0, -mdl.mass_T*mdl.CMx, 0;...
     0, -mdl.mass_b*mdl.CMz_b, 0, mdl.Ixx_b, 0, -mdl.Ixz_b;...
     mdl.mass_b*mdl.CMz_b, 0, -mdl.mass_T*mdl.CMx, 0, mdl.Iyy, 0;...
     0, mdl.mass_T*mdl.CMx, 0, -mdl.Ixz_b, 0, mdl.Izz]);

%********************Centers of Drag******************************
%calculate the centers of drag for the main components about the CG 
%start with the assumption that the center of drag is the center of
%the element
mdl.CDx_rotor = mdl.Xcg(25);
mdl.CDy_rotor = mdl.Ycg(25);
mdl.CDz_rotor = mdl.Zcg(25);
mdl.CDx_VBC = mdl.Xcg(26); 
mdl.CDy_VBC = 0;
mdl.CDz_VBC = mdl.Zcg(26);
mdl.CDx_body = (-(164.625/2-19.230)*0.0254-x_origin)*(20/3);
mdl.CDy_body = (0)*0.0254*(20/3);
mdl.CDz_body = (0)*0.0254*(20/3);

%*****************Important Locations, Distances, and Coefficients********

%Important dimensions for calculating drag
mdl.L_body = 146.605*0.0254*(20/3); %entire main body for drag calculations
mdl.D_body = 23.500*0.0254*(20/3);
mdl.L_VBC = 20; %set this
mdl.VBC_VB = 62.5;
VBC_total_disp_vol = mdl.VBC_VB + mdl.disp_mass(26)/mdl.rho; %total VBC displacement including entrained water
mdl.D_VBC = 2*sqrt(VBC_total_disp_vol/mdl.L_VBC/pi);

% Calculate the centers of the Cable attachment Point
mdl.CAPx = ((19.230)*0.0254 - x_origin)*(20/3)-5; %***********************
mdl.CAPy = 0*(20/3);
mdl.CAPz = mdl.D_body/2+1.011;
disp('Net positive buoyancy when tank is filled with air in N:')
disp((mdl.disp_mass_T+mdl.VBC_VB*mdl.rho-mdl.mass_T)*mdl.g)
disp('Net positive buoyancy when tank is filled with water in N:')
disp((mdl.disp_mass_T-mdl.mass_T)*mdl.g)

disp('Calculated length of main body in front of rotor in meters:')
disp(mdl.Xcg(5))
disp('Diameter of main body in meters:')
disp(mdl.D_body)
disp('Displaced Volume of VBC: Var Disp, Const Disp, Total Dist in m^3:')
disp(mdl.VBC_VB)
disp(mdl.disp_mass(26)/mdl.rho)
disp(VBC_total_disp_vol)
disp('Length of VBC in in meters:')
disp(mdl.L_VBC)
disp('Diameter of VBC in meters:')
disp(mdl.D_VBC)
disp('The CAP x, y, and z locations in meters are:')
disp(mdl.CAPx)
disp(mdl.CAPy)
disp(mdl.CAPz)

%Coeficients of drag
mdl.Cdx_body = 0.4;
mdl.Cdy_body = 1;
mdl.Cdz_body = 1;
mdl.Cdx_VBC = 0.2;
mdl.Cdy_VBC = 1;
mdl.Cdz_VBC = 1;

%cg buoyancy compensation modules (BCMs)
mdl.Xcg_VBC = mdl.Xcg(26);
mdl.Ycg_VBC = mdl.Ycg(26);
mdl.Zcg_VBC = mdl.Zcg(26);

