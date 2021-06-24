function waves = waves_spec_calc2(n,m,Hs,angle,s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAVE_SPEC_CALC2                     %
                                 %
%   % n is the number of discrete sampling intervals for w 
%   % m is the number of discrete sampling intervals for teta 
%   % Hs is the significant wave height     
%   % angle is the main direction of the spectrum
%   % s is the spreading parameter of the directional spectrum                                 %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Number_of_frequencies=n;
wmin=0.2;
wmax=2;

Number_of_directions=m;
tetamin=pi/180*(-90);
tetamax=pi/180*(90);

g=9.81;
h=320;
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

teta1=pi/180*angle;

w_step=(wmax-wmin)/(Number_of_frequencies-1);
teta_step=(tetamax-tetamin)/(Number_of_directions-1);

w=wmin:w_step:wmax; %wave frequencies
teta=tetamin:teta_step:tetamax;

k=w.^2/g; 

Number_of_waves=Number_of_frequencies*Number_of_directions;

pri=rand(1,Number_of_waves)*2*pi; %random phase given to each wave

%%%%%%%%Spectrum about the direction
N=1/pi*2^(2*s-1)*(factorial(s))^(2)/(factorial(2*s));
D=N.*(cos((teta-teta1)./2)).^(2*s);
%%%%%%%%

%%%%%%%%Spectrum about the frequency
    alpha=0.0081*g^2;
    beta=3.11/(Hs^2);
    S=alpha*w.^(-5).*exp(-beta.*w.^(-4));   %spectrum
%%%%%%%%

[S_prim,D_prim]=meshgrid(S,D);
A=sqrt(2.*S_prim.*w_step.*D_prim.*teta_step);%Wave amplitude of each wave

H=[];
w2=[];
teta2=[];
k2=[];

for i=1:Number_of_directions
    H=[H,A(i,:)];
    w2=[w2 w];
    k2=[k2 k];
    teta2=[teta2 ones(1,Number_of_frequencies).*teta(i)];
end

waves.teta2=teta2;
waves.g=g;
waves.h=h;
waves.H=H;
waves.pri=pri;
waves.w2=w2;
waves.k2=k2;

% %%%%%%%%%%%%%%
% %%%%%%%%%%%%%%