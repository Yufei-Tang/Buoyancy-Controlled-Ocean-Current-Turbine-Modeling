function [out] = calc_cable_water_vel(in,Mod_vars,ADCP)
Z=in(1:end-1);
t=in(end);

ADCP.inflowSelection = 1;
%water total water velocity in the inertial coordinate system
%InflowSelection=1 if using an imported dataset else use uniform values.
if ADCP.inflowSelection == 1 
    time_in = ADCP.tStart+datenum([0,0,0,0,0,t]);
    Uw_interp = interp2(ADCP.z,ADCP.t,ADCP.u',Z,time_in);
    Vw_interp = interp2(ADCP.z,ADCP.t,ADCP.v',Z,time_in);
    Ww_interp = interp2(ADCP.z,ADCP.t,ADCP.w',Z,time_in);       
else  
    % This is the water velicity if the use defines the surface water
    % velocity and the shear
    Uw_interp = Mod_vars.Current(1)+Mod_vars.Current(4)*Z;
    Vw_interp = Mod_vars.Current(2)+Mod_vars.Current(5)*Z;
    Ww_interp = Mod_vars.Current(3)+Mod_vars.Current(6)*Z;
end

out = [Uw_interp Vw_interp Ww_interp];



