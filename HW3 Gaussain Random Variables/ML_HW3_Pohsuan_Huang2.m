% Machine Learnign Exercise 3
% Problem 2(d)
% Po-Hsuan Huang 2014.11.14
% Plot vairance of Xh as well as posterior variance as a function of
% variance of auditory data. 

close all;  % close all figures
%clc;       % clear command area
clear;     % clear variables


%%  plot the variances
sigma= 4.9:0.001:5.1;   %figure 2b2
%sigma= 4:10;            %figure 2b3
%sigma= 5:50;            %figure 2b1


meanv = -10;
varv = 5^2;
meana = 10;
vara = sigma.^2;
mean0 = 0; % mean of prior
var0 = 100^2; % variance of prior

alpha =vara.^-1;
beta =varv^-1;
gamma = var0^-1;


% Xh =Xa+Xv
varh = (1/4)*(vara+varv);
%posterior mean
varpost = (alpha+beta+gamma).^-1;
% min(vara+varv)
minvar = min(vara,varv);



figure(1)
plot(sigma, varpost,sigma, varh,sigma,minvar);
title(  'posterior variance and variance of Xh');
legend('poster variance','varianc of Xh','min-variance' );
xlabel('sigma of auditory data');
ylabel('value');




%%when will Variance of h be worse than variance of visual signal?

   
   
        
    
   


