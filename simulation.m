clear 
clc
Amn = 1;    % Particular mode Constant
% A10 = 1;    % for example
a=input('Enter the width of the waveguide(m):'); %0.0229
b=a/2;
z=40*10^-3;
M = 40; %number of steps    
% % f = 450 *10^9;
% m=1;
% n=0;
f_temp = input('Enter operating ferquency(Ghz) f:');
f=f_temp*10^9;
cond = input('Enter the conductivity of the walls:');
%cond=5.8*10^7 ;                  % the conductivity for copper
epsilono = 8.8540e-12;           % Permittivity constant
epsilon_r = input('Enter the relative permeativity of the medium:');                  % Relative Permittivity constant
eps = epsilono*epsilon_r;
muo = (pi)*4e-7;               % Permeability constant
mu1_r = input('Enter the relative permeability of the medium:'); 
mu = muo*mu1_r;
omega = 2*pi*f; 
k=omega*sqrt(mu*eps);
Rs=sqrt((omega*muo)/(2*cond));
% kx=m*pi/a;
% ky=n*pi/b;
% kc=sqrt(kx^2+ky^2);
% beta=sqrt(k^2-kc^2);
% fc = kc/2*pi*sqrt(eps*mu);
mode=input('Enter choice: 1 for TE and 2 for TM:    ');

XScale = linspace(0,a,M);
YScale = linspace(0,b,M);
ZScale = linspace(0,z,M);


[x, y ,z]=meshgrid(XScale, YScale, ZScale);

%TE

if mode==1
m = input('Enter mode value m:');
n = input('Enter mode value n:');
% m=1;
% n=0;
kx=m*pi/a;
ky=n*pi/b;
kc=sqrt(kx^2+ky^2);
disp(kc);
beta=sqrt(k^2-kc^2);
fc = kc/(2*pi*sqrt(eps*muo));

disp(fc);
if f<=fc
    disp('The operating frequency is less than the cutoff frequency try again.');
else
Ex = ((1i*omega*mu*n*pi)/(kc^2*b))*Amn*cos(kx.*x).*sin(ky.*y).*exp(-1i*beta*z);

Ey = -((1i*omega*mu*m*pi)/(kc^2*a))*Amn*sin(kx.*x).*cos(ky.*y).*exp(-1i*beta*z);

Ez = 0;

Ez = zeros(size(real(Ey))); %??

Hx = ((1i*beta*m*pi)/(kc^2*a))*Amn*sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-1i*beta*z);

Hy = ((1i*n*pi)/(kc^2*b))*Amn*cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-1i*beta*z);

Hz = Amn*cos(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-1i*beta*z);
figure()
subplot(1,2,1);
quiver3 (x, y, z, real(Ex), real(Ey), real(Ez));
subplot(1,2,2);
quiver3 (x, y, z, real(Hx), real(Hy), real(Hz),'r');
figure()
hold on
quiver3 (x, y, z, real(Ex), real(Ey), real(Ez));
quiver3 (x, y, z, real(Hx), real(Hy), real(Hz),'r');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('E-H Graph');


if m==1&&n==0 
F=linspace(0,2*fc,M);
Omega = 2*pi*F;
K=Omega*sqrt(eps*mu);
RRS=sqrt((Omega*muo)/(2*cond));
Kc=ones(1,M)*kc;
Beta=sqrt(K.^2-Kc.^2);
PT=(mu*a^3*Amn^2*b*real(Beta).*Omega)/(4*pi^2);
PL=Amn^2*(b*ones(1,M)+(a/2)*ones(1,M)+Beta.^2*a^3/(2*pi^2) ).*RRS;
Attenuation=PL./(2*PT);

figure
plot(F,Attenuation);
xlabel('Frequency (Hz)');
ylabel('Attenuation due to conductor loss (Np/m)');
title('The attenuation with respect to Freqency in TE10 mode');
end

end
%TM MODE%
elseif mode ==2
Bmn=1;
m = input('Enter mode value m:');
n = input('Enter mode value n:');

% m=5;
% n=5;
kx=m*pi/a;
ky=n*pi/b;
kc=sqrt(kx^2+ky^2);
beta=sqrt(k^2-kc^2);
fc = kc/(2*pi*sqrt(eps*mu));

if f<=fc
    disp('The operating frequency is less than the cutoff frequency try again.');
else


Ex = ((-1i*beta*m*pi)/(kc^2*a))*Amn*cos(kx.*x).*sin(ky.*y).*exp(-1i*beta*z);

Ey = -((1i*beta*n*pi)/(kc^2*b))*Amn*sin(kx.*x).*cos(ky.*y).*exp(-1i*beta*z);

Ez = Bmn*sin(kx.*x).*cos(ky.*y).*exp(-1i*beta*z);


Hx = ((1i*omega*eps*n*pi)/(kc^2*b))*Bmn*sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-1i*beta*z);

Hy = ((-1i*omega*m*eps*pi)/(kc^2*a))*Bmn*cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-1i*beta*z);

Hz = zeros(size(real(Hy)));

figure()
hold on
quiver3 (x, y, z, real(Ex), real(Ey), real(Ez));
quiver3 (x, y, z, real(Hx), real(Hy), real(Hz),'r');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('E-H Graph');
 figure()
 subplot(1,2,1);
 quiver3 (x, y, z, real(Ex), real(Ey), real(Ez));
  subplot(1,2,2);
 quiver3 (x, y, z, real(Hx), real(Hy), real(Hz),'r');
end

end