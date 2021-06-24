function [out] = equations_of_motion_cable(in,Mod_vars)

out=[];
N_cable=Mod_vars.N_cable;
Mass_node=2*Mod_vars.m/N_cable; %2 for added mass

for i=1:(N_cable-1)
X_dot=in(i*6-5);
Y_dot=in(i*6-3);
Z_dot=in(i*6-1);

FX=in((N_cable-1)*6+i*3-2);
FY=in((N_cable-1)*6+i*3-1);
FZ=in((N_cable-1)*6+i*3);

Xv_dot = FX/Mass_node;
Yv_dot = FY/Mass_node;
Zv_dot = FZ/Mass_node;

out = [out;Xv_dot;X_dot;Yv_dot;Y_dot;Zv_dot;Z_dot];
end




