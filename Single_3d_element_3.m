function [sys,x0,str,ts] = Single_3d_element_3(t,x,u,flag,Mod_vars,ADCP)
%Mod_vars.Z0_ZCAP=depth of deepest node
N_cable=Mod_vars.N_cable;

if length(u) == (N_cable+1)*6+1
    if u((N_cable+1)*6+1) == 1
        x0 = 1;
        str = 1;
        ts = 1;
    end
end

switch flag

    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0
        [sys,x0,str,ts]=mdlInitializeSizes(Mod_vars,N_cable,ADCP);

        %%%%%%%%%%%%%%%
        % Derivatives %
        %%%%%%%%%%%%%%%
    case 1
        sys=mdlDerivatives(t,x,u,Mod_vars);

        %%%%%%%%%%
        % Update %
        %%%%%%%%%%
    case 2
        sys=mdlUpdate(t,x,u);

        %%%%%%%%%%%
        % Outputs %
        %%%%%%%%%%%
    case 3
        sys=mdlOutputs(t,x,u,Mod_vars,N_cable,ADCP);

        %%%%%%%%%%%%%%%%%%%%%%%
        % GetTimeOfNextVarHit %
        %%%%%%%%%%%%%%%%%%%%%%%
    case 4
        sys=mdlGetTimeOfNextVarHit(t,x,u);

        %%%%%%%%%%%%%
        % Terminate %
        %%%%%%%%%%%%%
    case 9
        sys=mdlTerminate(t,x,u);

        %%%%%%%%%%%%%%%%%%%%
        % Unexpected flags %
        %%%%%%%%%%%%%%%%%%%%
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);

end


function [sys,x0,str,ts]=mdlInitializeSizes(Mod_vars,N_cable,ADCP)

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = (N_cable+1)*3;
sizes.NumInputs      = (N_cable+1)*6+1;
sizes.DirFeedthrough = 6;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)

%
sys = [];

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u,Mod_vars,N_cable,ADCP)



% input variables defined
% u(1) velocity of node 1 (fixed one) in the x direction in the inertial frame
% u(2) position of node 1 in the x direction in the inertial frame
% u(3) velocity of node 1 in the y direction in the inertial frame
% u(4) position of node 1 in the y direction in the inertial frame
% u(5) velocity of node 1 in the z direction in the inertial frame
% u(6) position of node 1 in the z direction in the inertial frame

% u(i*6-5) velocity of node i (intermediate one) in the x direction in the inertial frame
% u(i*6-4) position of node i in the x direction in the inertial frame
% u(i*6-3) velocity of node i in the y direction in the inertial frame
% u(i*6-2) position of node i in the y direction in the inertial frame
% u(i*6-1) velocity of node i in the z direction in the inertial frame
% u(i*6) position of node i in the z direction in the inertial frame

% u((N_cable+1)*6-5) velocity of node N_cable+1 (Turbine one) in the x direction in the inertial frame
% u((N_cable+1)*6-4) position of node N_cable+1 in the x direction in the inertial frame
% u((N_cable+1)*6-3) velocity of node N_cable+1 in the y direction in the inertial frame
% u((N_cable+1)*6-2) position of node N_cable+1 in the y direction in the inertial frame
% u((N_cable+1)*6-1) velocity of node N_cable+1 in the z direction in the inertial frame
% u((N_cable+1)*6) position of node N_cable+1 in the z direction in the inertial frame ***


sys=zeros((N_cable+1)*3,1);

for i=1:N_cable %for each section of the cable

    x1=u(i*6-4);
    y1=u(i*6-2);
    z1=u(i*6);
    
    x2=u((i+1)*6-4);
    y2=u((i+1)*6-2);
    z2=u((i+1)*6);
    
    vx1=u(i*6-5);
    vy1=u(i*6-3);
    vz1=u(i*6-1);
    
    vx2=u((i+1)*6-5);
    vy2=u((i+1)*6-3);
    vz2=u((i+1)*6-1);

    f = zeros(6,1);

    L = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);

    %%  calculating the rotation matrix for element i
    if (z2-z1 == 0)
        teta = pi/2*sign(x2-x1);
    else
        teta = atan2((x2-x1),(z2-z1));
    end

    if L == 0
        fi = 0;
    else
        fi = asin(-(y2-y1)/L);
    end

    LIB(1,1) = cos(teta);
    LIB(1,2) = sin(teta)*sin(fi);
    LIB(1,3) = sin(teta)*cos(fi);
    LIB(2,1) = 0;
    LIB(2,2) = cos(fi);
    LIB(2,3) = -sin(fi);
    LIB(3,1) = -sin(teta);
    LIB(3,2) = cos(teta)*sin(fi);
    LIB(3,3) = cos(fi)*cos(teta);

    LBI = LIB';
    % transfering velocity of the nodes into the body fixed frame
    
    %%ADCP addition
    InflowSelection = 1;
    %water total water velocity in the inertial coordinate system
    %InflowSelection=1 if using an imported dataset else use uniform values.
    if InflowSelection == 1 
        time_in = ADCP.tStart+datenum([0,0,0,0,0,t]);
        z = [z1 z2];
        Uw_interp_x1 = interp2(ADCP.z,ADCP.t,ADCP.u',z1,time_in);
        Vw_interp_v1 = interp2(ADCP.z,ADCP.t,ADCP.v',z1,time_in);
        Ww_interp_w1 = interp2(ADCP.z,ADCP.t,ADCP.w',z1,time_in);       
        Uw_interp_x2 = interp2(ADCP.z,ADCP.t,ADCP.u',z2,time_in);
        Vw_interp_v2 = interp2(ADCP.z,ADCP.t,ADCP.v',z2,time_in);
        Ww_interp_w2 = interp2(ADCP.z,ADCP.t,ADCP.w',z2,time_in);
        vel1 = LBI*([vx1-Uw_interp_x1 vy1-Vw_interp_v1 vz1-Ww_interp_w1]'); 
        vel2 = LBI*([vx2-Uw_interp_x2 vy2-Vw_interp_v2 vz2-Ww_interp_w2]');
    else  
        % This has been modified to account for water velocity on 7/12/2020
        vel1 = LBI*([vx1-Mod_vars.Current(1)-Mod_vars.Current(4)*z1 vy1-Mod_vars.Current(2)-Mod_vars.Current(5)*z1 vz1-Mod_vars.Current(3)-Mod_vars.Current(6)*z1]');
        vel2 = LBI*([vx2-Mod_vars.Current(1)-Mod_vars.Current(4)*z2 vy2-Mod_vars.Current(2)-Mod_vars.Current(5)*z2 vz2-Mod_vars.Current(3)-Mod_vars.Current(6)*z2]');
    end

    % Internal tension
    if ((L - Mod_vars.l/N_cable) > 0)
        f(3) = Mod_vars.K*(L-Mod_vars.l/N_cable);
    else
        f(3) = 0;
    end

    f(6) = -f(3);

    %% Internal Damping
    f(3) = f(3) + Mod_vars.Cid/Mod_vars.l*(vy2 - vy1);
    f(6) = f(6) - Mod_vars.Cid/Mod_vars.l*(vx2 - vx1);

    % calculating hydrodynamic loading functions to be used in the drag
    % calculations

    [fn1,ft1,fn2,ft2] = cal_hydro_load(vel1,vel2);

    %total speed in the body fixed frame, upper node
    v2 = norm(vel2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% DRAG %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%% Upper node %%%%%%%%%%%%%%%%

    Velinnorplane2 = norm(vel2(1:2));

    if (abs(Velinnorplane2) >= eps)
        RatioB0 = vel2(1)/Velinnorplane2;
        RatioB1 = vel2(2)/Velinnorplane2;
        DragB0 = -Mod_vars.drag_n*fn2*(v2*v2)*(RatioB0);
        DragB1 = -Mod_vars.drag_n*fn2*(v2*v2)*(RatioB1);
    else
        DragB0 = 0;
        DragB1 = 0;
    end

    f(4) = f(4) + DragB0;
    f(5) = f(5) + DragB1;
      
    DragB2 = -Mod_vars.drag_t*ft2*v2*v2*sign(vel2(3));

    f(6) = f(6) + DragB2;
       
    sys((i*3-2):(i*3),1)= sys((i*3-2):(i*3),1) + [LIB*f(1:3,1)];
    sys(((i+1)*3-2):((i+1)*3),1)= sys(((i+1)*3-2):((i+1)*3),1) + [LIB*f(4:6,1)];
 
    %%%%%%%Weight of cable, accounting for the sea floor%%%%%%%%%%%
    if (i+1)<(N_cable+1)
        if u(i*6)<=u(6) %this is true if the cable node being evaluated is higher than the anchor
            sys((i+1)*3)=sys((i+1)*3)+Mod_vars.Wt/N_cable;
        elseif u(i*6)<=u(6)+4 %this is true if the cable node being evaluated is below the anchor by less than 1 m.
            sys((i+1)*3)=sys((i+1)*3)+Mod_vars.Wt/N_cable*(1-(u(i*6)-u(6))/2);
        else %this if for the nod being more than 1 m below the anchor where it is modeled as being nutrally buoyant
            sys((i+1)*3)=sys((i+1)*3)-Mod_vars.Wt/N_cable;
        end      
    else
        sys((i+1)*3)=sys((i+1)*3)+Mod_vars.Wt/N_cable/2;
    end
   
end

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate


%************************************************************
%************************************************************

% calculating the tension all cables except the first
function f = Internal_Forces(i,x,Mod_vars)
% i is the element number

L = x(2*(i+1))-x(2*i);   % length of first element
dl = L-Mod_vars.l; % strech of first element

% check to see if it is in tension
if dl > 0
    % if the cable is in tension, then we calculated the tension due to
    % elastic stretch and internal damping
    f = dl*Mod_vars.K+Mod_vars.Cid*(x(2*(i+1)-1)-x(2*i-1));
else
    % if the cable is slack, we set the tension to zero
    f = 0;
end


%************************************************************
%************************************************************

% calculating the drag force on the on the cable nodes
function f = Drag_forces(i,x,Mod_vars)

% i is the node number
V = 2/3*x(2*i-1)*abs(x(2*i-1))+1/6*x(2*i-3)*abs(x(2*i-3))+1/6*x(2*i+1)*abs(x(2*i+1)) ...
    -1/12*(x(2*i-3)-x(2*i-1))*abs(x(2*i-3)-x(2*i-1))-1/12*(x(2*i-1)-x(2*i+1))*abs(x(2*i-1)-x(2*i+1));

%f = 805.0331*V;
f = Mod_vars.drag_c*V;


function [fn1,ft1,fn2,ft2] = cal_hydro_load(vr1,vr2)

if (vr1(3) == 0)
    beta1 = 1.57079632679;
else
    beta1 = atan(sqrt(vr1(1)*vr1(1)+vr1(2)*vr1(2))/vr1(3));
end
if (vr2(3) == 0)
    beta2 = 1.57079632679;
else
    beta2 = abs(atan(sqrt(vr2(1)*vr2(1)+vr2(2)*vr2(2))/vr2(3)));
end

%fn1 and ft1 are normal and tangential loading functions at the upper node
fn1 = (0.5 + (-0.1)*cos(beta1) + 0.1*sin(beta1) + (-0.4)*cos(2*beta1) + (-.011)*sin(2*beta1));
ft1 = (-0.1945 + 0.203*cos(beta1) + 0.1945*sin(beta1) + 0*cos(2*beta1)...
    + (-0.0681)*sin(2*beta1)-(-0.0681)*sin(pi));

%fn2 and ft2 are normal and tangential loading functions at the lower node
fn2 = (0.5 + (-0.1)*cos(beta2) + 0.1*sin(beta2) + (-0.4)*cos(2*beta2) + (-.011)*sin(2*beta2));
ft2 = (-0.1945 + 0.203*cos(beta2) + 0.1945*sin(beta2) + 0*cos(2*beta2) +...
    (-0.0681)*sin(2*beta2)-(-0.0681)*sin(pi));
