function [lny_m,lny_sd] = attenuation(M,R,V)

%Based on the paper:
%Kenneth W. Campbell and Yousef Bozorgnia, "NGA-West2 Ground Motion Model for the 
%Average Horizontal Components of PGA, PGV, and 5% Damped Linear Acceleration Response Spectra", 
%Earthquake Spectra, Volume 30, No. 3, pages 1087–1115, 2014.

%only for M>6.5; V>865m/s

%model parameters of NGA West-2
c0= -4.416 ; c1= .984; c2=0.537 ;
c3=-1.499 ; c4=-.496 ; c5=-2.773 ; c6=.248 ;
c7 =6.768 ; c11=1.090; k1=865 ; k2=-1.186 ; n=1.18 ;

f_mag = c0+c1*M+c2*(M-4.5)+c3*(M-5.5)+c4*(M-6.5); %magnitude term
f_dist = (c5+c6*M)*log(sqrt(R^2+c7^2)); %distance from source term
f_soil = (c11 + k2*n)*log(V/k1); %soil conditions term
%Note: original NGA West-2 model has more terms that are considered 0 here

lny_m = f_mag+f_dist+f_soil; % lnPGA - mean value of normal
lny_sd = sqrt(.322^2+.492^2); % standard deviation of lnPGA random variable
