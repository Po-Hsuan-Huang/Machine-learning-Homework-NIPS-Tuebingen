% Machine Learnign Exercise 3
% Problem 2(b)
% Po-Hsuan Huang 2014.11.14
% Plot posterior as well as the (visual and auditory)likehood funcitons
%for mean0 = 0, variance0=100, meanv =-10,meana= 10, variancev = 5, variancea= 20. 

close all;  % close all figures
%clc;       % clear command area
clear;     % clear variables

%% plot the visual likelihood function

x= -40:40;

meanv = -10;
varv = 5^2;

Visual = (sqrt(2*pi*varv)^-1) * exp(-((x-meanv).^2)/(2*varv));



%% plot the auditory likelihood function

meana = 10;
vara = 20^2;

Auditory = (sqrt(2*pi*vara)^-1) * exp(-((x-meana).^2)/(2*vara));

%% plot the posterior probability function

mean0 = 0; % mean of prior
var0 = 100^2; % variance of prior

alpha =vara^-1;
beta =varv^-1;
gamma = var0^-1;


meanpost = (alpha*meana+beta*meanv+gamma*mean0)/(alpha+beta+gamma);  

varpost = (alpha+beta+gamma)^-1;

Probability =@(x) sqrt(2*pi*varpost).^-1*exp(-((x-meanpost).^2)/2) ; 
% joint P(x)= likelihood*prior

Normal= integral(Probability, -Inf,Inf);  %Normalization constant for posterior

Postdist = sqrt(2*pi*varpost).^-1*exp(-((x-meanpost).^2)/2)/Normal;
%% plot all 
figure(1)
subplot(3,1,1)
plot(x,Visual);
title(  'visual likelihood function');
xlabel('x');
ylabel('p(x)');
subplot(3,1,2)
plot(x,Auditory);
title(  'auditory likelihood function');
xlabel('x');
ylabel('p(x)');
subplot(3,1,3)
plot(x,Postdist);
title(  'Posterior distribution');
xlabel('x');
ylabel('p(x)');


figure(2)
plot(x,Visual);
hold on
plot(x,Auditory);
plot(x,Postdist);
title(  'postdistribution and its likelihood functions');
legend('visual','autitoy','posterior' );
xlabel('x');
ylabel('p(x)');
hold off

