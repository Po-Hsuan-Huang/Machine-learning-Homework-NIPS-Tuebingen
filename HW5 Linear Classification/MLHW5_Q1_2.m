% Machine Learning Exercise 5
% Problem 2(a) ,(b)
% Po-Hsuan Huang 2014.11.26
% calculate covariances and means for the given xTrain.
% classify 2 classes of Gaussian distribution with linear discriminant.




function MLHW5_Q1_2(tTrain,xTrain,tTest,xTest)
%% ALways clean the mess first
close all;  % close all figures
% control +C kills the process


%% linear discriminant



% for t=1 group
C1X= xTrain(tTrain==1,1);
C1Y= xTrain(tTrain==1,2);
mean1= [mean(C1X) mean(C1Y)];
cov1 = cov(C1X, C1Y);   % calculate the covariant matrix.

det(cov1);
% for t=-1 group
C2X= xTrain(tTrain==-1,1);
C2Y= xTrain(tTrain==-1,2);
mean2= [mean(C2X) mean(C2Y)];
cov2 = cov(C2X,C2Y);   % calculate the covariant matrix
det(cov2);
cov_avg= (cov1+cov2)/2;   % mean of the two covariant matrix

%calculate weight and threshold.

Wg=  -(mean2-mean1)/cov_avg; 

W0 = -(mean1*inv(cov_avg)*mean1'- mean2*inv(cov_avg)*mean2')/2;

display(Wg);
display(W0);
%% mis classifying rate

% trainig set
misCL = MCR(C1X,C1Y,C2X,C2Y,Wg,W0);
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

cov_avg_D= (cov1_D+cov2_D)/2;   % mean of the two covariant matrix


%% you dont need this part becuase you should use the test set to 
% evaluate the proformance of your discrimanant garnered from training set.

% 
% Wg_D=  -(mean2_D-mean1_D)/cov_avg_D; 
% 
% W0_D = -(mean1_D*inv(cov_avg_D)*mean1_D'- mean2_D*inv(cov_avg_D)*mean2_D')/2;

misCL_D = MCR(D1X,D1Y,D2X,D2Y,Wg,W0);
display(misCL_D);
%% plotting  Traning set
figure(1)
plot(C1X,C1Y,'ro',C2X,C1Y,'bx');


title('Q 2-B linear classifying');
xlabel('xTrain(:,1)');
ylabel('xTrain(:,2)');
hold on
%  classifing plane =     Wg*s'+omega0 = 0;
x= -10:10;
y =   -(Wg(1)*x+W0)/Wg(2);
plot(x,y); 
legend('Class1','Class2')

hold off

%% plotting Test set

figure(2)
plot(D1X,D1Y,'ro',D2X,D2Y,'bx');
title('Q 2-B linear classifying')
xlabel('xTest(:,1)');
ylabel('xTest(:,2)');
hold on
%  classifing plane =     Wg*s'+omega0 = 0;
x= -10:10;
y =   -(Wg_D(1)*x+W0_D)/Wg_D(2);
plot(x,y); 
legend('Class1','Class2')


hold off


end

function [misCL]=MCR(C1X,C1Y,C2X,C2Y,Wg,W0)
  %% misclassifying rate
     C1(1,:)=C1X; 
     C1(2,:)=C1Y;
     C2(1,:)=C2X;
     C2(2,:)=C2Y;
    for  count = 1 :size(C1X)
        
        C1_Err(count) =   Wg*C1(:,count)+W0<=0;
        C2_Err(count) =   Wg*C2(:,count)+W0>=0;
    end      
  
% for class1 
misCL_1=  mean(C1_Err);
% for class2
misCL_2=   mean(C2_Err);
%  for total

    
misCL = (sum(C1_Err)+sum(C2_Err))/500;





end