% Machine Learning Exercise 3
% Problem 2(2)
% Po-Hsuan Huang 2014.11.20
% Training Data with linear regression model.

%% ALways clean the mess first
close all;  % close all figures
% control +C kills the process

%% use alpha = beta= 1, calculate the posterior mean

% basically this is a totally irrelavent questions from above.
sigma= 5;
feature= zeros(51,20);
for i = 1:50
    feature(1,:)= xTrain';
    feature(i+1,:)=  2*exp(-(xTrain-i).^2/sigma^2);
    
end

%% Calculate the Nx50 matrix ztrain for which the n-th row is  Z(Xn) and produce 
%the image of the matrix.

zTrain = feature(2:51,:)';


figure(2)
%From the for loop generating the following plot we can tell 
%each vertical stripe corresponding to the correlation between training 
%set xTrain with basis function i

imagesc(zTrain);
ylabel('xTrain data points');
 xlabel('basis function');

alpha = 1;
beta = 1; 
Covpost= inv(alpha*eye(50)+beta*zTrain'*zTrain);  % covariance matrix
Meanpost= beta*Covpost*zTrain'*tTrain;             % posterior mean


figure(3)
plot( Meanpost,'-*');
title('Q2_b')
ylabel('PostMean');
xlabel('basis function');



%% Q2_3 weighted predictive mean using xspace 0:0.2:50
  feature2= zeros(51,251);
  xspace = 0:0.2:50;
  feature2(1,:)= xspace; 
  Awesum= zeros(1,251);
  figure(5)

 for i=1:50   
    feature2(i+1,:)=  Meanpost(i).*2*exp(-(xspace-i).^2/sigma^2);
    Awesum(1,:)= Awesum(1,:) +  Meanpost(i).*2*exp(-(xspace-i).^2/sigma^2);
    plot(xspace,feature2(i+1,:));
    
hold on
 end   
 plot(xTrain,tTrain,'*');
 plot(xspace,Awesum(1,:));
title('Q2_c')
ylabel('t');
xlabel('x');
legend('predictive mean', 'tTrain');
hold off

figure()
surf(feature2(2:51,:),'EdgeColor','none');
xlabel('basis function')
ylabel('x space')
title('basis function weighting w influnce on predictive mean at x space ')
colorbar

%% Predict mean using xTrain
feature3= zeros(51,20);
  
  feature3(1,:)= xTrain; 
  Meanpredict= zeros(1,20);
  figure(5)

 for i=1:50   
    feature3(i+1,:)=  Meanpost(i).*2*exp(-(xTrain'-i).^2/sigma^2);
    Meanpredict(1,:)= Meanpredict(1,:) +  Meanpost(i).*2*exp(-(xTrain'-i).^2/sigma^2);
    plot(xTrain,feature3(i+1,:));
    
hold on
 end   
 plot(xTrain,tTrain,'*');
 plot(xTrain,Meanpredict(1,:));
title('Q2_c')
ylabel('t');
xlabel('x');
legend('predictive mean', 'tTrain');
hold off


%% Q2_4 display covariance as image 
figure(6)
surf(Covpost,'EdgeColor','none');
title('posterior covariance metrix,Q2_d');
sddv_post=sqrt(diag(Covpost));
figure(7)
errorbar(Meanpost,2*sddv_post,'-o');
axis('tight')
title('Q2_d_2')
ylabel('PostMean');
xlabel('basis function');

%% Predictive variance is sum(Zi*Varpost)
% for i = 1:50
%     
%     feature2(3,:)= Varpost(i)*2*exp(-(xspace-i).^2/sigma^2);
% end
% Varpredict= feature2(3,:);
% sddv_predict= sqrt(Varpredict);
% figure(8)
% errorbar(Meanpredic,2*sddv_predict);
% title('Q2_e')
% ylabel('x');
% xlabel('t');
Varpost = diag(Covpost);

for i =1:20  
   Varpredict(1,:) = 1/beta+ Varpost(i)*(exp(-(xTrain-i).^2/sigma^2)).^2;
end

%Varpredict= feature2(3,:);
sddv_predict= sqrt(Varpredict);
figure(8)
errorbar(xTrain,Meanpredict,2*sddv_predict,'-o');
title('Q2_e')
ylabel('t');
xlabel('x');
hold on
plot(xTrain,tTrain,'*');
legend('predictive mean', 'tTrain');
hold off
