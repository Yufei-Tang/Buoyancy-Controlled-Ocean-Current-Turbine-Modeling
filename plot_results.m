clear all
clc

load 7dofTI5
torqueTI5=Torque;
omegaTI5=Omega;
powerTI5=246480*omegaTI5(:,2)*2*pi/60;
positionTI5=position;
ang_velTI5=ang_velocity;
eulerTI5=euler;
turb_forcesTI5=turb_forces;
fx_onebladeTI5= fx_oneblade;
clear Torque Omega position ang_velocity euler turb_forces fx_oneblade

load 7dofTI5C5fpoint001
torqueTI5_f0=Torque;
omegaTI5_f0=Omega;
powerTI5_f0=246480*omegaTI5_f0(:,2)*2*pi/60;
positionTI5_f0=position;
ang_velTI5_f0=ang_velocity;
eulerTI5_f0=euler;
turb_forcesTI5_f0=turb_forces;
fx_onebladeTI5_f0= fx_oneblade;
clear Torque Omega position ang_velocity euler turb_forces fx_oneblade

load 7dofTI5C27
torqueTI5_c27=Torque;
omegaTI5_c27=Omega;
powerTI5_c27=246480*omegaTI5_c27(:,2)*2*pi/60;
positionTI5_c27=position;
ang_velTI5_c27=ang_velocity;
eulerTI5_c27=euler;
turb_forcesTI5_c27=turb_forces;
fx_onebladeTI5_c27= fx_oneblade;
clear Torque Omega position ang_velocity euler turb_forces fx_oneblade

load 7dofTI20
torqueTI20=Torque;
omegaTI20=Omega;
powerTI20=246480*omegaTI20(:,2)*2*pi/60;
positionTI20=position;
ang_velTI20=ang_velocity;
eulerTI20=euler;
turb_forcesTI20=turb_forces;
fx_onebladeTI20= fx_oneblade;
clear Torque Omega position ang_velocity euler turb_forces fx_oneblade

figure (1)
subplot(1,2,1)
plot(torqueTI5(:,1)-60,torqueTI5(:,2)/1000,'r',torqueTI5(:,1)-60,powerTI5(:,1)/1000,':b')
legend('torque','power')
xlabel('time (s)'); ylabel('shaft power (kW)/hydrodynamic torque (kN-m)')
axis([0 max(torqueTI5(:,1))-60 140 475])
title('TI=5%')
subplot(1,2,2)
plot(torqueTI20(:,1)-60,torqueTI20(:,2)/1000,'r',torqueTI20(:,1)-60,powerTI20(:,1)/1000,':b')
xlabel('time (s)'); 
axis([0 max(torqueTI5(:,1))-60 140 475])
title('TI=20%')

avg_power5=mean(powerTI5(1201:end,1))
dev_power5=std(powerTI5(1201:end,1))
avg_power20=mean(powerTI20(1201:end,1))
dev_power20=std(powerTI20(1201:end,1))

avg_torque5=mean(torqueTI5(1201:end,2))
dev_torque5=std(torqueTI5(1201:end,2))
avg_torque20=mean(torqueTI20(1201:end,2))
dev_torque20=std(torqueTI20(1201:end,2))

figure (2)
subplot(1,2,1)
plot(fx_onebladeTI5(:,1)-60,-fx_onebladeTI5(:,2)/1000,'b')
title('TI=5%')
xlabel('time (s)'); ylabel('Axial force on one blade(kN)')
axis([0 max(torqueTI5(:,1))-60 100 190])
subplot(1,2,2)
plot(fx_onebladeTI20(:,1)-60,-fx_onebladeTI20(:,2)/1000,'b')
title('TI=20%')
xlabel('time (s)'); 
axis([0 max(torqueTI5(:,1))-60 100 190])

avg_axialForce5=mean(fx_onebladeTI5(1201:end,2))
dev_axialForce5=std(fx_onebladeTI5(1201:end,2))
avg_axialForce20=mean(fx_onebladeTI20(1201:end,2))
dev_axialForce20=std(fx_onebladeTI20(1201:end,2))

figure (3)
subplot(1,2,1)
plot(eulerTI5(:,1)-60,eulerTI5(:,2)*180/pi,'r',eulerTI5(:,1)-60,eulerTI5(:,3)*180/pi,'--g',eulerTI5(:,1)-60,eulerTI5(:,4)*180/pi,':b')
legend('roll','pitch','yaw')
title ('TI=5%')
xlabel('time (s)'); ylabel('degrees')
axis([0 max(torqueTI5(:,1))-60 -3 1])
subplot(1,2,2)
plot(eulerTI20(:,1)-60,eulerTI20(:,2)*180/pi,'r',eulerTI20(:,1)-60,eulerTI20(:,3)*180/pi,'--g',eulerTI20(:,1)-60,eulerTI20(:,4)*180/pi,':b')
title ('TI=20%')
xlabel('time (s)'); 
axis([0 max(torqueTI5(:,1))-60 -3 1])

figure (4)
subplot(1,2,1)
plot(fx_onebladeTI5(:,1)-60,-fx_onebladeTI5(:,2)/1000,'k',fx_onebladeTI5_c27(:,1)-60,-fx_onebladeTI5_c27(:,2)/1000,'--r')
legend('C=5','C=27')
xlabel('time (s)'); ylabel('Axial force on one blade(kN)')
axis([0 max(torqueTI5(:,1))-60 110 150])
subplot(1,2,2)
plot(torqueTI5(:,1)-60,powerTI5(:,1)/1000,'k',torqueTI5(:,1)-60,powerTI5_c27(:,1)/1000,'--r')
xlabel('time (s)'); ylabel('Shaft power(kW)')
axis([0 max(torqueTI5(:,1))-60 330 410])

figure (5)
subplot(1,2,1)
plot(fx_onebladeTI5(:,1)-60,-fx_onebladeTI5(:,2)/1000,'k',fx_onebladeTI5_f0(:,1)-60,-fx_onebladeTI5_f0(:,2)/1000,'--r')
legend('f_m_i_n=0.01','f_m_i_n=0.001')
xlabel('time (s)'); ylabel('Axial force on one blade(kN)')
axis([0 max(torqueTI5(:,1))-60 100 160])
subplot(1,2,2)
plot(torqueTI5(:,1)-60,powerTI5(:,1)/1000,'k',torqueTI5_f0(:,1)-60,powerTI5_f0(:,1)/1000,'--r')
xlabel('time (s)'); ylabel('Shaft power(kW)')
axis([0 max(torqueTI5(:,1))-60 330 410])

figure (6)
plot(positionTI5(1201:end,3),positionTI5(1201:end,4),':b',positionTI20(1201:end,3),positionTI20(1201:end,4),'--r')
legend('TI=5%','TI=20%')
xlabel('Y location (m)'); ylabel('Z location (m)')

figure (7)
subplot (1,3,1)
plot(positionTI20(1201:end,1),positionTI20(1201:end,2))
xlabel('time (s)'); ylabel('X location (m)')
subplot (1,3,2)
plot(positionTI20(1201:end,1),positionTI20(1201:end,3))
xlabel('time (s)'); ylabel('Y location (m)')
subplot (1,3,3)
plot(positionTI20(1201:end,1),positionTI20(1201:end,4))
xlabel('time (s)'); ylabel('Z location (m)')
% 
% figure (6)
% subplot (1,2,1)
% plot(position_ti0(1201:end,1),position_ti0(1201:end,3))
% xlabel('time (s)'); ylabel('Y location (m)')
% subplot (1,2,2)
% plot(position_ti0(1201:end,1),position_ti0(1201:end,4))
% xlabel('time (s)'); ylabel('Z location (m)')





