% Machine Learning Exercise 5
% Problem 2(e)
% Po-Hsuan Huang 2014.11.26
% Use guadratic discriminant ayalysis to classify 2 class.
% without assuming the two covariance are the same.
% plot the decision boundary on the scattering plot of data.



function MLHW5_Q1_2_C(tTrain,xTrain,tTest,xTest)
%% ALways clean the mess first
close all;  % close all figures
% control +C kills the process





%% Quadratic discriminant

% for t=1 group
C1X= xTrain(tTrain==1,1);  
C1Y= xTrain(tTrain==1,2);
% The logic judgement in the first arguement only maintain those elements
% in xTrain whose corresponding tTrain value is 1.
% Therefore, the length of C1X becomes 250.

mean1= [mean(C1X) mean(C1Y)];
cov1 = cov(C1X, C1Y);   % calculate the covariant matrix.

det1 = det(cov1);
% for t=-1 group
C2X= xTrain(tTrain==-1,1);
C2Y= xTrain(tTrain==-1,2);
mean2= [mean(C2X) mean(C2Y)];
cov2 = cov(C2X,C2Y);   % calculate the covariant matrix
det2= det(cov2);



%calculate A,B,C for discrminant xAx' + Bx' + C for training phase.

A_c=  (inv(cov2)-inv(cov1))/2; 

B_c = (mean1/(cov1)- mean2/(cov2));

C_c = (mean2/(cov2)*mean2' - mean1/(cov1)*mean1')/2 - log(det1/det2);

%% mis classifying rate

% trainig set
misCL = MCR(C1X,C1Y,C2X,C2Y,A_c,B_c,C_c);
display(misCL);
%% test set


% for t=1 group
D1X= xTest(tTest==1,1);
D1Y= xTest(tTest==1,2);
mean1_D= [mean(D1X) mean(D1Y)];
cov1_D =cov(D1X,D1Y);   % calculate the covariant matrix.


% for t=-1 group
D2X= xTest(tTest==-1,1);
D2Y= xTest(tTest==-1,2);
mean2_D= [mean(D2X) mean(D2Y)];
cov2_D = cov(D2X,D2Y);   % calculate the covariant matrix


% The following part should only exert to training phase, not test phase!!
% Or the xTest data would be another training data instead of testing data.

%calculate A,B,C for discrminant xAx' + Bx' + C for training phase.
% 
% A_d=  (inv(cov2_D)-inv(cov1_D))/2; 
% 
% B_d = (mean1_D/(cov1_D)- mean2_D/(cov2_D));
% 
% C_d = (mean2_D/(cov2_D)*mean2_D' - mean1_D/(cov1_D)*mean1_D')/2;



misCL_D = MCR(D1X,D1Y,D2X,D2Y,A_c,B_c,C_c);
display(misCL_D);



%% plotting  Traning set
figure(1)
plot(C1X,C1Y,'ro',C2X,C2Y,'bx');
hold on
 
title('Q 2-E quadratic classifying')
xlabel('xTrain(:,1)');
ylabel('xTrain(:,2)');

%classifing plane =     x'Ax + b'x + c  = 0;
[x1, x2 ] = meshgrid(-10:10,-10:10);

% It's hard to do it matric wise,
%z= [x1,x2]'*A_c*[x1,x2] +B_c*[x1,x2]+C_c;

% do it element wise, sigh...
z= A_c(1,1)*x1.^2+(A_c(1,2)+A_c(2,1)).*x1.*x2+A_c(2,2).*x2.^2 + B_c(1).*x1+B_c(2).*x2+C_c;

contour(x1,x2,z,'ShowText','on')  
legend('Class1','Class2')
hold off


%% plotting Test set

figure(2)
plot(D1X,D1Y,'ro',D2X,D2Y,'bx');
hold on
title('Q 2-E quadratic classifying')
xlabel('xTest(:,1)');
ylabel('xTest(:,2)');


%classifing plane =     x'Ax + b'x + c  = 0;
[x1, x2 ] = meshgrid(-10:10,-10:10);

% It's hard to do it matric wise,
%z= [x1,x2]'*A_c*[x1,x2] +B_c*[x1,x2]+C_c;

% do it element wise, sigh...
z= A_c(1,1)*x1.^2+(A_c(1,2)+A_c(2,1)).*x1.*x2+A_c(2,2).*x2.^2 + B_c(1).*x1+B_c(2).*x2+C_c;

contour(x1,x2,z,'ShowText','on')  
legend('Class1','Class2')
hold off
 


end

function [misCL]=MCR(C1X,C1Y,C2X,C2Y,A,B,C)
  %% misclassifying rate
     C1(1,:)=C1X; 
     C1(2,:)=C1Y;
     C2(1,:)=C2X;
     C2(2,:)=C2Y;
    for  count = 1 :size(C1X)
        
        C1_Err(count) =   (C1(:,count)'*A*C1(:,count)+B*C1(:,count)+C)<=0;
        C2_Err(count) =   (C2(:,count)'*A*C2(:,count)+B*C2(:,count)+C)>=0;
    end    

% misCL_1=  mean(C1_Err);
% % for class2
% misCL_2=   mean(C2_Err);
% %  for total
% 
% misCL = (misCL_1+misCL_2)/2; 



misCL = sum(C1_Err+C2_Err)/(length(C1)+length(C2));





end