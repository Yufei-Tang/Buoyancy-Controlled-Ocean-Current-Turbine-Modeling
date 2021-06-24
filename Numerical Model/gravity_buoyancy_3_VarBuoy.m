function OUT = gravity_buoyancy_3_VarBuoy(states, mdl)

%u = states(1);%not used
%v = states(2);%not used 
%w = states(3);%not used 
%p_b = states(4);%not used
%p_r = states(5);%not used
%q = states(6);%not used 
%r = states(7);%not used 
%X = states(8);%not used
%Y = states(9);%not used
%Z = states(10);%not used
phi_b = states(11);   %roll body
%phi_r = states(12);   %roll rotor
theta = states(13); %pitch
psi = states(14);   %yaw
tank_fill_front = states(15);  %front tank fill fraction (values between 0 and 1) ***
tank_fill_back = states(16);  %back tank fill fraction (values between 0 and 1) ***

phi = phi_b;
%Transforms from body to inertial coordinates
L_IB = [[cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)];...
        [cos(theta)*sin(psi) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi)];...
        [-sin(theta)         sin(phi)*cos(theta)                            cos(phi)*cos(theta)                           ]]; 

%Transforms from inertial to body coordinates
L_BI = L_IB';

%%%%%%%%%%%%%% GRAVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%forces in the inertial frame
f_XYZ_G_b = [0 0 mdl.mass_b*mdl.g]';
f_XYZ_G_r = [0 0 mdl.mass_r*mdl.g]';
%calculate the forces in the body fixed frame
f_xyz_G_b = L_BI*f_XYZ_G_b;
f_xyz_G_r = L_BI*f_XYZ_G_r;
f_xyz_G = f_xyz_G_r+f_xyz_G_b;
M_xyz_G_b = cross([mdl.CMx_b-mdl.CMx mdl.CMy_b-mdl.CMy mdl.CMz_b-mdl.CMz],f_xyz_G_b);
M_xyz_G_r = cross([mdl.CMx_r-mdl.CMx mdl.CMy_r-mdl.CMy mdl.CMz_r-mdl.CMz],f_xyz_G_r);%moment about center of gravity
%M_xyz_G = M_xyz_G_b+M_xyz_G_r;     %moment about center of gravity

%%%%%%%%%%%%%% BUOYANCY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%center of buoyancy for front tank in the x direction = mdl.Xcg_BCM + 0.25*mdl.L_ear
%center of buoyancy for back tank in the x direction = mdl.Xcg_BCM - 0.25*mdl.L_ear
%center of buoyancy for both tanks in the y direction = mdl.Ycg_BCM
%center of buoyancy for both tanks in the z direction = mdl.Zcg_BCM
%the entire original displaced mass of both BSMs combined = mdl.disp_mass(26)
disp_mass_VarBuoy_front = mdl.disp_mass(26)*tank_fill_front; %***
disp_mass_VarBuoy_back = mdl.disp_mass(26)*tank_fill_back; %***

%forces in the inertial frame
f_XYZ_B_VarBuoy_front = [0 0 -disp_mass_VarBuoy_front*mdl.g]'; %***
f_XYZ_B_VarBuoy_back = [0 0 -disp_mass_VarBuoy_back*mdl.g]'; %***
f_XYZ_B_b = [0 0 -mdl.disp_mass_b*mdl.g]';
f_XYZ_B_r = [0 0 -mdl.disp_mass_r*mdl.g]';

%calculate the forces in the body fixed frame
f_xyz_B_VarBuoy_front = L_BI*f_XYZ_B_VarBuoy_front; %***
f_xyz_B_VarBuoy_back = L_BI*f_XYZ_B_VarBuoy_back; %***
f_xyz_B_b = L_BI*f_XYZ_B_b;
f_xyz_B_r = L_BI*f_XYZ_B_r;
f_xyz_B = f_xyz_B_b + f_xyz_B_r + f_xyz_B_VarBuoy_front + f_xyz_B_VarBuoy_back; %***


M_xyz_B_VarBuoy_front = cross([(mdl.Xcg_VBC + 0.25*mdl.L_VBC)-mdl.CMx mdl.Ycg_VBC-mdl.CMy mdl.Zcg_VBC-mdl.CMz],f_xyz_B_VarBuoy_front);     %moment about center of gravity
M_xyz_B_VarBuoy_back = cross([(mdl.Xcg_VBC - 0.25*mdl.L_VBC)-mdl.CMx mdl.Ycg_VBC-mdl.CMy mdl.Zcg_VBC-mdl.CMz],f_xyz_B_VarBuoy_back);     %moment about center of gravity
M_xyz_B_b = cross([mdl.CBx_b-mdl.CMx mdl.CBy_b-mdl.CMy mdl.CBz_b-mdl.CMz],f_xyz_B_b);     %moment about center of gravity
M_xyz_B_r = cross([mdl.CBx_r-mdl.CMx mdl.CBy_r-mdl.CMy mdl.CBz_r-mdl.CMz],f_xyz_B_r);     %moment about center of gravity

f_xyz = f_xyz_G + f_xyz_B;
M_x_r = M_xyz_G_r(1) + M_xyz_B_r(1);
M_x_b = M_xyz_G_b(1) + M_xyz_B_b(1) + M_xyz_B_VarBuoy_front(1) + M_xyz_B_VarBuoy_back(1);
M_yz = M_xyz_G_r(2:3) + M_xyz_G_b(2:3) + M_xyz_B_r(2:3) + M_xyz_B_b(2:3) + M_xyz_B_VarBuoy_front(2:3) + M_xyz_B_VarBuoy_back(2:3);


OUT = [f_xyz; M_x_b; M_x_r; M_yz'];