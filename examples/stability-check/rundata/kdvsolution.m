%%%%%%%%%%%%%%%%%%
% kdvsolution.m
% calculate ekdv and kdv solution for internal waves
% yun zhang @stanford
%%%%%%%%%%%%%%%%%%%
% parameter
h1=10;
h2=990;
rho1=999.5;
rho2=1000.5;
g=9.81;
eta0=0.1;
% prepared variable
sigma=2*(rho2-rho1)/(rho2+rho1);
c0=sqrt(g*sigma*h1*h2/(h1+h2));
alpha1=3/2*c0*(h2-h1)/(h1*h2);
alpha2=3*c0/(h1^2*h2^2)*(7/8*(h1-h2)^2-(h2^3+h1^3)/(h1+h2));
alpha2=0;
beta=c0*h1*h2/6;
% wave property
c=c0+eta0/3*(alpha1+0.5*alpha2*eta0);
L=1/sqrt(eta0*(alpha1+0.5*alpha2*eta0)/12/beta);
b=-eta0*alpha2/(2*alpha1+alpha2*eta0);
L %843.27

