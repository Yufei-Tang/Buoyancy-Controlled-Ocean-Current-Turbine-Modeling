function turbulence = turbulence_constants_v3(TI,fmin,fmax,radius,N_mesh,U,a,b,C,P)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% turbulence_constants_v2                    %
                                 %
%   % input(1): TI is the turbulence intencity 
%   % input(2): minimum frequency, fmin
%   % input (3): maximum frequency, fmax                            %                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Put all of your math here
dth=2*pi/N_mesh; %azimuthal angle
theta=linspace(dth,2*pi,N_mesh);

%Calculation of components of TI
sigma=U*TI;% Overall rms of all fluctuating components, sigma=sqrt((sigma_x^2+sigma_y^2+sigma_z^2))
sigma_x=sqrt(1/(1+a^2+b^2))*sigma; %rms of fluctuating component in x-direction
sigma_y=a*sigma_x; %rms of fluctuating component in y-direction
sigma_z=b*sigma_x; %rms of fluctuating component in z-direction
TIx=sigma_x/U; %TI in x-direction
TIy=sigma_y/U; %TI in y-direction
TIz=sigma_z/U; %TI in z-direction, TI=sqrt((TIx^2+TIy^2+TIz^2)/3)

for i=1:length(radius)
   zp=-radius*cos(theta);
   yp=radius*sin(theta);
end

Pz=(zp(:))';
Py=(yp(:))';

%this one is for the rotor center variable buoyancy chamber center
Pz=[0 Pz -15.34];
Py=[0 Py 0];

M=length(Pz); %Points in space

% distance among points
for i=1:M
    for j=1:M
        rx=(Pz(i)-Pz(j));
        ry=(Py(i)-Py(j));
        r(i,j)=sqrt(rx^2+ry^2);
    end
end


%turbulence simulation
df=(fmax-fmin)/P;
f=fmin:df:fmax; %frequency discretized
Ax=4/3*TIx^2*(U^2)/((1/(fmin)^(2/3)-1/(fmax)^(2/3))); %proportionality constant in x-direction
Ay=4/3*TIy^2*(U^2)/((1/(fmin)^(2/3)-1/(fmax)^(2/3))); %proportionality constant in y-direction
Az=4/3*TIz^2*(U^2)/((1/(fmin)^(2/3)-1/(fmax)^(2/3))); %proportionality constant in z-direction

for i=1:P
    fm(i)=(f(i)+f(i+1))/2;
    Coh=exp(-C*r*fm(i)/U);
    Sx=Coh*Ax*fm(i)^(-5/3)*df; %PSD in x-direction
    Sy=Coh*Ay*fm(i)^(-5/3)*df; %PSD in y-direction
    Sz=Coh*Az*fm(i)^(-5/3)*df; %PSD in z-direction
    Hx=chol(Sx, 'lower');%Cholesky Decomposition in x-direction
    Hy=chol(Sy, 'lower');%Cholesky Decomposition in y-direction
    Hz=chol(Sz, 'lower');%Cholesky Decomposition in z-direction
    ph1=(exp(1i*repmat(2*pi*rand(M,1),1,M)'));
    ph2=(exp(1i*repmat(2*pi*rand(M,1),1,M)'));
    ph3=(exp(1i*repmat(2*pi*rand(M,1),1,M)'));
    hx=Hx.*ph1;
    hy=Hy.*ph2;
    hz=Hz.*ph3;
    vel_inter_x = sum(hx,2);
    vel_inter_y = sum(hy,2);
    vel_inter_z = sum(hz,2);
    for j=1:M
        vel_f_x(i,j)=vel_inter_x(j);
        vel_f_y(i,j)=vel_inter_y(j);
        vel_f_z(i,j)=vel_inter_z(j);
        mag_x(i,j) = abs(vel_f_x(i,j));
        mag_y(i,j) = abs(vel_f_y(i,j));
        mag_z(i,j) = abs(vel_f_z(i,j));
        ang_x(i,j) = 2*pi+angle(vel_f_x(i,j));
        ang_y(i,j) = 2*pi+angle(vel_f_y(i,j));
        ang_z(i,j) = 2*pi+angle(vel_f_z(i,j));
    end
end

r_uv=-0.136*sigma_x;%cross correlation between u&v, user input
r_uw=-0.079*sigma_x-0.325;%cross correlation between u&w, user input

%% This is where you decide what you will be passing to your primary code
turbulence.ang_x=ang_x;
turbulence.ang_y=ang_y;
turbulence.ang_z=ang_z;
turbulence.mag_x=mag_x;
turbulence.mag_y=mag_y;
turbulence.mag_z=mag_z;
turbulence.fm=fm;
turbulence.r_uv=r_uv;
turbulence.r_uw=r_uw;
turbulence.M=M;
turbulence.P=P;
turbulence.U=U;