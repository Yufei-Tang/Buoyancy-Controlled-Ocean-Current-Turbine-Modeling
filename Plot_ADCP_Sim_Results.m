clc
clear
%close all
load('TestADCP_09_25_2017_at_1800hr_v1.mat')
%% User Inputs Go Here
%This is the minutes of data that 
min_start = 10;
tStart = datenum('25-Sep-0017 18:00:00');%start time EST for generating a subset of data v1
tStop = datenum('27-Sep-0017 12:00:00');%stop time EST for generating a subset of data
%% Simulation Calculations Go Here
start_date = tStart;
%This should calculate the datenum for all data
date_num = start_date+datenum(0,0,0,0,0,out.euler(:,1));
%This should calculate the number of data points removed from the start of
%the numerical simulation
dt = out.euler(2,1)-out.euler(1,1);
num_start = round(60*min_start/dt);
%This should calculate the power used by the ballast system
d_FF1 = [0; diff(out.FF(:,2))]; %this is the change in FF each time step
d_FF2 = [0; diff(out.FF(:,3))];
neg_FF1 = d_FF1.*floor(d_FF1); %this is the removal of ballast water only
neg_FF2 = d_FF2.*floor(d_FF2);
E_empty = 0.01402; %Amount of energy required to completely empty a ballast tank at a depth of 50 m in MW*hr. 
P_ballast = E_empty*3600/dt*(neg_FF1+neg_FF2);
%% ADCP Calculations Go Here
filename = 'ADCP-Data-Irma.mat'; %Input file with ADCP data
ADCP = ADCP_Init(filename,tStart,tStop,1);

%% Plotting Code Goes Here
figure(1)
%This subplot should show the current speed
subplot(5,1,1)
pcolor(ADCP.t-5/24,-ADCP.z,(ADCP.u.^2+ADCP.v.^2).^0.5);
caxis([0 2.5]);colorbar('vert');shading('flat');
axis([date_num(num_start)-5/24 date_num(end)-5/24 -300 0])
title(strcat('Start Date & Time: ',datestr(start_date-5/24,21),' GMT'))
ylabel('Current Speed [m/s]')
datetick('x',15,'keeplimits');

%This subplot should show the current direction
subplot(5,1,2)
pcolor(ADCP.t-5/24,-ADCP.z,180/pi*atan2(ADCP.v,ADCP.u));
caxis([-40 40]);colorbar('vert');shading('flat');
axis([date_num(num_start)-5/24 date_num(end)-5/24 -300 0])
ylabel('Current Direction [deg]')
datetick('x',15,'keeplimits')

%This subplot should show the Euler angles
subplot(5,1,3)
plot(date_num(num_start:end)-5/24,out.euler(num_start:end,2:3)*180/pi,date_num(num_start:end)-5/24,(pi-out.euler(num_start:end,4))*180/pi)
axis tight
ylabel('Euler Angles [deg]')
legend('\phi','\theta','\psi','Location','best')
datetick('x',15,'keeplimits')

%This subplot should show the NED positions 
subplot(5,1,4)
plot(date_num(num_start:end)-5/24,out.position(num_start:end,2:4))
axis tight
ylabel('Location - NED [m]')
legend('North','East','Down','Location','best')
datetick('x',15,'keeplimits')

%This should plot the produced and consumed powers
subplot(5,1,5)
plot(date_num(num_start:end)-5/24,-out.tau_EM(num_start:end,2).*out.Omega(num_start:end,2)/1000000,date_num(num_start:end)-5/24,P_ballast(num_start:end))
axis tight
xlabel('Time [HH:MM]')
ylabel('power [MW]')
legend('Produced','Ballast')
datetick('x',15,'keeplimits')

figure(11)
%This subplot should show the current speed
pcolor(ADCP.t-5/24,-ADCP.z,(ADCP.u.^2+ADCP.v.^2).^0.5);
caxis([0 2.5]);colorbar('vert');shading('flat');
axis([date_num(num_start)-5/24 date_num(end)-5/24 -300 0])
title(strcat('Start Date & Time: ',datestr(start_date-5/24,21),' GMT'))
xlabel('Time [HH:MM]')
ylabel('Current Speed [m/s]')
datetick('x',15,'keeplimits');

figure(12)
%This subplot should show the current direction
pcolor(ADCP.t-5/24,-ADCP.z,180/pi*atan2(ADCP.v,ADCP.u));
caxis([-40 40]);colorbar('vert');shading('flat');
axis([date_num(num_start)-5/24 date_num(end)-5/24 -300 0])
title(strcat('Start Date & Time: ',datestr(start_date-5/24,21),' GMT'))
xlabel('Time [HH:MM]')
ylabel('Current Direction [deg]')
datetick('x',15,'keeplimits')

figure(13)
%This subplot should show the Euler angles
plot(date_num(num_start:end)-5/24,out.euler(num_start:end,2:3)*180/pi,date_num(num_start:end)-5/24,(pi-out.euler(num_start:end,4))*180/pi)
axis tight
title(strcat('Start Date & Time: ',datestr(start_date-5/24,21),' GMT'))
xlabel('Time [HH:MM]')
ylabel('Euler Angles [deg]')
legend('\phi','\theta','\psi','Location','best')
datetick('x',15,'keeplimits')

figure(14)
%This subplot should show the NED positions 
plot(date_num(num_start:end)-5/24,out.position(num_start:end,2:4))
axis tight
title(strcat('Start Date & Time: ',datestr(start_date-5/24,21),' GMT'))
xlabel('Time [HH:MM]')
ylabel('Location - NED [m]')
legend('North','East','Down','Location','best')
datetick('x',15,'keeplimits')

figure(15)
%This should plot the produced and consumed powers
plot(date_num(num_start:end)-5/24,-out.tau_EM(num_start:end,2).*out.Omega(num_start:end,2)/1000000,date_num(num_start:end)-5/24,P_ballast(num_start:end))
axis tight
title(strcat('Start Date & Time: ',datestr(start_date-5/24,21),' GMT'))
xlabel('Time [HH:MM]')
ylabel('power [MW]')
legend('Produced','Ballast')
datetick('x',15,'keeplimits')
%text(1.5e4,-2,'power [MW]')

