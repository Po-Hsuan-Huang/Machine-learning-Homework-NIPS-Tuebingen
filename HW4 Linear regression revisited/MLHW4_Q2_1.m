% Machine Learning Exercise 3
% Problem 2(1)
% Po-Hsuan Huang 2014.11.20
% Training Data with linear regression model.
%% ALways clean the mess first
close all;  % close all figures
% control +C kills the process

%% Construct feature vector with basis function 2exp(-(x-i)^2/sigma^2)
sigma= 5;
feature= zeros(51,501);
for i = 1:50
    feature(1,:)= xPlot';
    feature(i+1,:)=  2*exp(-(xPlot-i).^2/sigma^2);
    figure(1)
    plot(feature(1,:),feature(i+1,:),'.');
    hold on
end

title('Gaussian basis functions');
xlabel('x');
ylabel('f(x)');
hold off

%% Calculate the 50xN matrix ztrain for which the n-th row is  Z(Xn) and produce 
%the image of the matrix.

zTrain = feature(2:51,:)';
figure(2)
clims= [0 2];
imagesc(zTrain',clims);
title('design matrix');
ylabel('Basis functions f');
xlabel('f(xn)');


